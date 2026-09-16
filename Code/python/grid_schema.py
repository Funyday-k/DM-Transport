"""Read version 1 of the minimal C++ phase-space grid metadata schema."""

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Mapping, Tuple, Union


_UNITS = {"r": "cm", "v": "cm/s", "mu": "1", "cell_phase_measure": "cm^6/s^3"}
_ESCAPE_THRESHOLD = {
    "energy_definition": "0.5*(v_cm_s^2-v_escape_cm_s(r_cm)^2)",
    "bound_condition": "E<0",
    "unbound_condition": "E>0",
    "zero_surface": "zero_phase_measure",
    "touching_cell": "classify_by_interior",
    "point_deposition": "containing_fv_cell_no_fractional_split",
}
_REQUIRED_KEYS = frozenset({
    "schema_version", "r_faces_cm", "v_faces_cm_s", "mu_faces",
    "flatten_order", "units", "quadrature_order_per_axis", "escape_threshold",
})


def _faces(value: Any, name: str, nonnegative: bool) -> Tuple[float, ...]:
    """Validate strictly increasing finite faces in cm, cm/s, or dimensionless mu."""
    if not isinstance(value, list) or len(value) < 2:
        raise ValueError(f"{name} must contain at least two faces")
    faces = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ValueError(f"{name} must contain JSON numbers")
        try:
            number = float(item)
        except OverflowError as error:
            raise ValueError(f"{name} must contain finite faces") from error
        if not math.isfinite(number) or (nonnegative and number < 0.0):
            raise ValueError(f"{name} must contain finite valid faces")
        if faces and number <= faces[-1]:
            raise ValueError(f"{name} must strictly increase")
        faces.append(number)
    return tuple(faces)


def _reject_constant(value: str) -> None:
    """Reject non-JSON NaN/Infinity tokens in metadata numeric values."""
    raise ValueError(f"nonfinite JSON constant: {value}")


def _unique_object(pairs: list) -> dict:
    """Reject duplicate JSON keys rather than silently changing a convention."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


@dataclass(frozen=True)
class GridSchema:
    """Versioned r [cm], v [cm/s], mu grid geometry; mu is fastest."""

    r_faces_cm: Tuple[float, ...]
    v_faces_cm_s: Tuple[float, ...]
    mu_faces: Tuple[float, ...]
    quadrature_order_per_axis: int

    @classmethod
    def from_mapping(cls, metadata: Mapping[str, Any]) -> "GridSchema":
        """Validate one complete version-1 metadata object with inline faces."""
        if not isinstance(metadata, Mapping) or set(metadata) != _REQUIRED_KEYS:
            raise ValueError("grid schema keys do not match version 1")
        if type(metadata["schema_version"]) is not int or metadata["schema_version"] != 1:
            raise ValueError("unsupported grid schema version")
        if metadata["flatten_order"] != ["r", "v", "mu"]:
            raise ValueError("flatten order must be r, v, mu from slowest to fastest")
        if metadata["units"] != _UNITS:
            raise ValueError("grid units do not match version 1")
        if metadata["escape_threshold"] != _ESCAPE_THRESHOLD:
            raise ValueError("escape-threshold convention does not match version 1")
        order = metadata["quadrature_order_per_axis"]
        if type(order) is not int or order not in (1, 2):
            raise ValueError("quadrature order per axis must be 1 or 2")
        r = _faces(metadata["r_faces_cm"], "r_faces_cm", True)
        v = _faces(metadata["v_faces_cm_s"], "v_faces_cm_s", True)
        mu = _faces(metadata["mu_faces"], "mu_faces", False)
        if mu[0] != -1.0 or mu[-1] != 1.0:
            raise ValueError("mu faces must span exactly [-1, 1]")
        return cls(r, v, mu, order)

    @property
    def shape(self) -> Tuple[int, int, int]:
        """Return (r, v, mu) cell counts, dimensionless."""
        return (len(self.r_faces_cm) - 1,
                len(self.v_faces_cm_s) - 1,
                len(self.mu_faces) - 1)

    def flatten(self, ir: int, iv: int, imu: int) -> int:
        """Map zero-based (r, v, mu) cell indices to mu-fastest flat index."""
        nr, nv, nmu = self.shape
        if any(type(index) is not int for index in (ir, iv, imu)) or not (
            0 <= ir < nr and 0 <= iv < nv and 0 <= imu < nmu
        ):
            raise IndexError("phase-space cell index is out of range")
        return (ir * nv + iv) * nmu + imu

    def cell_phase_measure(self, ir: int, iv: int, imu: int) -> float:
        """Integrate 8*pi^2*r^2*v^2 dr dv dmu in cm^6/s^3."""
        self.flatten(ir, iv, imu)
        r_lo, r_hi = self.r_faces_cm[ir:ir + 2]
        v_lo, v_hi = self.v_faces_cm_s[iv:iv + 2]
        mu_lo, mu_hi = self.mu_faces[imu:imu + 2]
        factors = (
            8.0 * math.pi * math.pi / 9.0,
            r_hi - r_lo, r_hi, r_hi, 1.0 + r_lo / r_hi + (r_lo / r_hi) ** 2,
            v_hi - v_lo, v_hi, v_hi, 1.0 + v_lo / v_hi + (v_lo / v_hi) ** 2,
            mu_hi - mu_lo,
        )
        mantissa, exponent = 1.0, 0
        for factor in factors:
            part, power = math.frexp(factor)
            mantissa *= part
            mantissa, adjustment = math.frexp(mantissa)
            exponent += power + adjustment
        try:
            measure = math.ldexp(mantissa, exponent)
        except OverflowError as error:
            raise OverflowError("phase-space cell measure overflows binary64") from error
        if not math.isfinite(measure) or measure <= 0.0:
            raise OverflowError("phase-space cell measure is not representable")
        return measure


def load_grid_schema(path: Union[str, Path]) -> GridSchema:
    """Read a standalone UTF-8 grid schema JSON file with cgs faces."""
    with Path(path).open("r", encoding="utf-8") as source:
        metadata = json.load(source, parse_constant=_reject_constant,
                             object_pairs_hook=_unique_object)
    return GridSchema.from_mapping(metadata)
