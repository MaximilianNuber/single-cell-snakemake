# scpipe/types.py
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, Mapping, Literal

Kind = Literal["embedding", "neighbors", "umap", "clusters", "qc"]

@dataclass(frozen=True)
class Result:
    kind: Kind
    outputs: Mapping[str, Any] = field(default_factory=dict)  # small arrays/tables only
    state:   Mapping[str, Any] = field(default_factory=dict)  # params, schema/version
    metrics: Mapping[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

# # yourpkg/types.py
# from __future__ import annotations

# from dataclasses import dataclass, field, asdict
# from typing import Any, Dict, Mapping, Optional

# @dataclass(frozen=True)
# class Result:
#     """
#     A tiny, immutable bundle returned by analysis functions.

#     - outputs: arrays/frames to be written into AnnData (e.g. {"layer": arr, "obs.size_factors": sf})
#     - state:   learned parameters / fitted values (e.g. {"target_sum": 1e4})
#     - metrics: scalar diagnostics (e.g. {"cv_sf": 0.12})

#     Keep it small and JSON-ish so it can be logged/serialized easily.
#     """
#     outputs: Mapping[str, Any] = field(default_factory=dict)
#     state: Mapping[str, Any] = field(default_factory=dict)
#     metrics: Mapping[str, float] = field(default_factory=dict)

#     def to_dict(self) -> Dict[str, Any]:
#         """Convenience for logging/serialization."""
#         return asdict(self)
