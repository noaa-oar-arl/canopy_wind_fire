"""
Run canopy-app using the namelist and load nc data.
"""
from __future__ import annotations

import contextlib
import json
import os
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

import f90nml
import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).parent.absolute()
assert HERE.name == "python"
assert (HERE.parent / "README.md").is_file() and (
    HERE.parent / "LICENSE"
).is_file(), "repo base"
REPO = HERE.parent


def _load_default_config() -> f90nml.Namelist:
    with open(REPO / "input" / "namelist.canopy") as f:
        config = f90nml.read(f)

    def as_abs(p_str: str) -> str:
        p0 = Path(p_str)
        if not p0.is_absolute():
            p = REPO / p0
        else:
            p = p0
        return p.as_posix()

    # Make input paths absolute in default config
    for k, v in config["filenames"].items():
        if isinstance(v, list):
            config["filenames"][k] = [as_abs(v_) for v_ in v]
        else:
            config["filenames"][k] = as_abs(v)

    return config


DEFAULT_CONFIG = _load_default_config()


@contextlib.contextmanager
def out_and_back(p: Path, *, finally_: Callable | None = None):
    """Context manager: change working directory to path `p` but return to original cwd after."""
    cwd = os.getcwd()
    os.chdir(p)
    try:
        yield
    finally:
        if finally_ is not None:
            finally_()
        os.chdir(cwd)


DEFAULT_POINT_INPUT = pd.read_csv(
    REPO / "input" / "point_file_20220701.sfcf000.txt", index_col=False
)

_TXT_STEM_SUFFS = {
    "wind": "_canopy_wind",
    "waf": "_waf",
    "eddy": "_eddy",
    "phot": "_phot",
    "bio": "_bio",
}


def run(
    *,
    config: dict[str, Any] | None = None,
    case_dir: Path | None = None,
    cleanup: bool = True,
    verbose: bool = False,
) -> xr.Dataset:
    """Run canopy-app, returning xarray dataset.

    Parameters
    ----------
    config
        Namelist options as a dict.
        Defaults used if not provided.
        (Don't need to change input path, we set to absolute path automatically.)
        If using relative paths, specify relative to the working directory.
    case_dir
        Directory for the case.
        A temporary directory is used if not provided.
        Replicates a bit of the structure (``input``, ``output`` subdirectories).
    cleanup
        Clean up (remove) the case directory after (successful) run.
    """
    # Ready the output directory
    if case_dir is None:
        import tempfile

        temp_dir = tempfile.TemporaryDirectory(prefix="canopy-app_")
        case_dir = Path(temp_dir.name)
    else:
        case_dir.mkdir(parents=True, exist_ok=True)

    # Create config
    full_config: f90nml.Namelist = deepcopy(DEFAULT_CONFIG)
    user_config = config or {}
    for section_name, sub_config in user_config.items():
        if section_name not in full_config:
            raise ValueError(f"Unexpected namelist section: {section_name!r}")
        for k, v in sub_config.items():
            assert not isinstance(v, dict)
            if k not in full_config[section_name]:
                raise ValueError(f"Unexpected namelist option: '{section_name}.{k}'")
            if verbose:
                print(f"'{section_name}.{k}' {full_config[section_name][k]} -> {v!r}")
            full_config[section_name][k] = v
    output_dir = case_dir / "output"
    output_dir.mkdir(exist_ok=True)
    ofp_stem = output_dir / "out"
    full_config["filenames"]["file_out"] = ofp_stem.relative_to(case_dir).as_posix()

    # Check input file(s)
    input_files_setting = full_config["filenames"]["file_vars"]
    if isinstance(input_files_setting, list):
        ifps_to_check = [Path(s) for s in input_files_setting]
    else:
        ifps_to_check = [Path(input_files_setting)]
    nc_outs = []
    for ifp in ifps_to_check:
        if not ifp.is_file():
            raise ValueError(f"Input file {ifp.as_posix()!r} does not exist")
        if ifp.suffix in {
            ".nc",
            ".nc4",
            ".ncf",
        }:  # consistent with sr `canopy_check_input`
            nc_out = True
            assert full_config["userdefs"]["infmt_opt"] == 0
        elif ifp.suffix in {".txt"}:
            nc_out = False
            assert full_config["userdefs"]["infmt_opt"] == 1
        else:
            raise ValueError(f"Unexpected input file extension: {ifp.suffix}")
        nc_outs.append(nc_out)
    if not len(set(nc_outs)) == 1:
        raise ValueError(
            f"Expected all input files to be of the same type (nc or txt). "
            f"filenames.file_vars: {input_files_setting}."
        )
    nc_out = nc_outs[0]
    if isinstance(input_files_setting, list):
        abs_path_strs = []
        for s in input_files_setting:
            ifp = Path(s)
            if not ifp.is_absolute():
                abs_path_strs.append(ifp.resolve().as_posix())
            else:
                abs_path_strs.append(ifp.as_posix())
        full_config["filenames"]["file_vars"] = abs_path_strs
    else:
        ifp = Path(input_files_setting)
        if not ifp.is_absolute():
            full_config["filenames"]["file_vars"] = ifp.resolve().as_posix()

    # Write namelist
    if verbose:
        print("Full namelist config:")
        print(full_config)
    config_dir = case_dir / "input"
    config_dir.mkdir(exist_ok=True)
    full_config.write(config_dir / "namelist.canopy", force=True)

    # Run canopy-app
    exe = REPO / "canopy"
    if not exe.is_file():
        raise RuntimeError("compile canopy-app first")
    cmd = [exe.as_posix()]
    with out_and_back(case_dir):
        subprocess.run(cmd, check=True, capture_output=not verbose)

    # Load nc
    if nc_out:
        # Should be just one file, even if multiple output time steps
        patt = f"{ofp_stem.name}*.nc"
        cands = sorted(output_dir.glob(patt))
        if not cands:
            raise ValueError(
                f"No matches for pattern {patt!r} in directory {output_dir.as_posix()!r}. "
                f"Files present are: {[p.as_posix() for p in output_dir.glob('*')]}."
            )
        if len(cands) > 1:
            print("Taking the first nc file only.")
        ds0 = xr.open_dataset(cands[0], decode_times=True)
        ds = (
            ds0.rename_dims(grid_xt="x", grid_yt="y")
            .swap_dims(level="z")
            .assign(time=ds0.time)
            .set_coords(["time", "z", "lat", "lon"])
            .squeeze()
        )
    else:
        dfs: list[pd.DataFrame] = []
        ifcans = [
            k
            for k, v in full_config["userdefs"].items()
            if k.startswith("ifcan") and v is True
        ]
        if not ifcans:
            raise RuntimeError("None of the 'ifcan' switches are on, so no text output.")

        # Read txt outputs
        for ifcan in ifcans:
            which = ifcan[len("ifcan") :]
            stem_suff = _TXT_STEM_SUFFS.get(which)
            if stem_suff is None:
                print(
                    f"warning: skipping {which!r} ({ifcan}) output since stem suffix unknown."
                )
                continue
            # NOTE: Separate file for each time
            patt = f"{ofp_stem.name}_*{stem_suff}.txt"
            cands = sorted(output_dir.glob(patt))
            if not cands:
                raise ValueError(
                    f"No matches for pattern {patt!r} in directory {output_dir.as_posix()!r}. "
                    f"Files present are: {[p.as_posix() for p in output_dir.glob('*')]}."
                )
            if verbose:
                print(f"detected output files for {ifcan}:")
                print("\n".join(f"- {p.as_posix()}" for p in cands))
            dfs_ifcan = []
            for cand in cands:
                df_t = read_txt(cand)
                df_t["time"] = np.datetime64(df_t.attrs["time"], "ns")
                dfs_ifcan.append(df_t)
            df = pd.concat(dfs_ifcan, ignore_index=True)
            df.attrs.update(which=which)
            df.attrs.update(df_t.attrs)
            dfs.append(df)

        # Merge
        units: dict[str, str] = {}
        dss = []
        for df in dfs:
            if {"time", "lat", "lon", "height"}.issubset(df.columns):
                ds_ = df.set_index(["time", "height", "lat", "lon"]).to_xarray().squeeze()
            elif {"time", "lat", "lon"}.issubset(df.columns):
                ds_ = df.set_index(["time", "lat", "lon"]).to_xarray().squeeze()
            else:
                raise ValueError(
                    "Expected df to have columns 'time', 'lat', 'lon' [,'height']. "
                    f"Got: {sorted(df)}."
                )
            units.update(df.attrs["units"])
            for vn in ds_.data_vars:
                assert isinstance(vn, str)
                ds_[vn].attrs["units"] = units[vn]
                if vn != "lad":  # in multiple groups
                    ds_[vn].attrs["group"] = df.attrs["which"]
            ds_.attrs.update(
                {k: v for k, v in df.attrs.items() if k not in {"which", "units"}}
            )
            dss.append(ds_)
        ds = xr.merge(dss, combine_attrs="no_conflicts")
        ds = ds.rename(height="z")
        if {"lat", "lon"}.issubset(ds.dims):
            ds = ds.rename_dims(lat="y", lon="x")
        # NOTE: lat/lon are 1-D in our examples, though 2-D in the example nc output

    if "lad" in ds.data_vars:
        ds = ds.set_coords("lad")

    # Store namelist settings
    ds.attrs["nml"] = str(full_config)
    ds.attrs["nml_json"] = json.dumps(full_config.todict())

    # Clean up
    if cleanup:
        shutil.rmtree(case_dir)

    return ds


def read_txt(fp: Path) -> pd.DataFrame:
    """Read canopy-app txt output file."""
    import re

    # Parse header lines
    with open(fp) as f:
        for i, line in enumerate(f):
            if i == 0:
                pattern = r" *time stamp\: *([0-9\.\:\-]*)"
                m = re.match(pattern, line)
                if m is None:
                    raise ValueError(
                        f"Unexpected file format. Line {i} failed to match regex {pattern!r}."
                    )
                time_stamp = pd.Timestamp(m.group(1))
            elif i == 1:
                pattern = r" *reference height, h\: *([0-9\.]*) m"
                m = re.match(pattern, line)
                if m is None:
                    raise ValueError(
                        f"Unexpected file format. Line {i} failed to match regex {pattern!r}."
                    )
                href = float(m.group(1))
            elif i == 2:
                pattern = r" *number of model layers\: *([0-9]*)"
                m = re.match(pattern, line)
                if m is None:
                    raise ValueError(
                        f"Unexpected file format. Line {i} failed to match regex {pattern!r}."
                    )
                nlay = int(m.group(1))
            elif i == 3:
                # Column names (some with units)
                heads = re.split(r" {2,}", line.strip())
                names: list[str] = []
                units: dict[str, str | None] = {}
                for head in heads:
                    m = re.fullmatch(r"([a-zA-Z_ ]+) \(([^)]+)\)", head)
                    if m is None:
                        u = None
                        name = head
                    else:
                        name, u = m.groups()
                    clean_name = name.lower().replace(" ", "_")
                    names.append(clean_name)
                    units[clean_name] = u
            else:
                break
        else:
            raise ValueError(
                "Unexpected file format. Expected 4 header lines followed by data."
            )

    df = pd.read_csv(
        fp,
        index_col=False,
        skiprows=4,
        header=None,
        delimiter=r"\s+",
        dtype=np.float32,
    )
    if len(names) != len(df.columns):
        raise RuntimeError(
            f"Unexpected file format. Detected columns names {names} ({len(names)}) "
            f"are of a different number than the loaded dataframe ({len(df.columns)})."
        )
    df = df.replace(np.float32(-9.0e20), np.nan)  # fill value, defined in const mod
    df.columns = names  # type: ignore[assignment]
    df.attrs.update(href=href, nlay=nlay, units=units, time=time_stamp)

    return df


def run_config_sens(
    cases: list[dict[str, Any]],
    base_dir: Path | None = None,
    cleanup: bool = True,
    verbose: bool = False,
) -> xr.Dataset:
    """Do multiple runs with different namelist configuration,
    returning single merged dataset.

    Parameters
    ----------
    cases
        List of namelist configurations in dict format,
        each one suitable for passing to :func:`run`.
    base_dir
        Base directory for the cases (case directories are created under).
        A temporary directory is used if not provided.
    """
    from collections import defaultdict

    if base_dir is None:
        import tempfile

        temp_dir = tempfile.TemporaryDirectory(prefix="canopy-app_")
        base_dir = Path(temp_dir.name)

    dss = []
    case_vars = defaultdict(list)
    for i, case in enumerate(cases):
        assert not ({"modlays", "modres"} & set(case["userdefs"]))

        for k, v in case["userdefs"].items():
            case_vars[k].append(v)

        if verbose:
            print(f"Running case {i+1}/{len(cases)}")
        ds = run(
            config=case,
            case_dir=base_dir / f"case_{i:0{len(str(len(cases) - 1))}}",
            cleanup=False,
            verbose=verbose,
        )
        dss.append(ds)

    ds = xr.concat(dss, dim="case", join="exact", combine_attrs="drop_conflicts")

    for k, v in case_vars.items():
        ds[k] = ("case", v)

    # Clean up
    if cleanup:
        shutil.rmtree(base_dir)

    return ds


def _k_sec(k: str) -> str:
    if k.startswith("file_"):
        sec = "filenames"
    else:
        sec = "userdefs"
    return sec


def config_cases(*, product: bool = False, **kwargs) -> list[dict[str, Any]]:
    """Helper to generate config dicts.

    Supply any namelist parameter as a keyword argument,
    with a single option or list of options
    (those with lists must be the same length,
    unless using `product`, which generates all possible combinations).
    """
    from collections import defaultdict

    if not product:
        sings = {}
        mults = {}
        for k, v in kwargs.items():
            if not (
                np.isscalar(v)
                or (type(v) is list and np.isscalar(v[0]))
                or (type(v) is list and type(v[0]) is list and np.isscalar(v[0][0]))
            ):  # FIXME: could be more stringent
                raise ValueError(
                    f"Expected {k!r} value to be "
                    f"scalar, list[scalar], or list[list[scalar]], got: {type(v)}."
                )

            if k == "file_vars":
                # Only support single time step runs for now
                if type(v) is list:
                    assert type(v[0]) is str
                    mults[k] = v
                else:
                    assert type(v) is str
                    sings[k] = v
                continue

            if (np.isscalar(DEFAULT_CONFIG[_k_sec(k)][k]) and np.isscalar(v)) or (
                type(DEFAULT_CONFIG[_k_sec(k)][k]) is list
                and type(v) is list
                and np.isscalar(v[0])
            ):
                sings[k] = v
            else:
                mults[k] = v

        max_len = max(len(v) for v in mults.values())
        for k, v in mults.items():
            if len(v) != max_len:
                raise ValueError(
                    "All lists should be the same length. "
                    f"Max length is {max_len} but {k!r} has length {len(v)}."
                )

        cases = []
        for i in range(max_len):
            case: defaultdict[str, dict[str, Any]] = defaultdict(dict)
            for k, v in sings.items():
                case[_k_sec(k)][k] = v
            for k, v in mults.items():
                case[_k_sec(k)][k] = v[i]
            cases.append(case)

    else:
        import itertools

        kwargs_ = deepcopy(kwargs)
        for k, v in kwargs_.items():
            if not isinstance(v, list):
                # NOTE: This condition won't work if the namelist default is an array
                kwargs_[k] = [v]

        cases = []
        for vs in itertools.product(*kwargs_.values()):
            case = defaultdict(dict)
            for k, v in zip(kwargs_, vs):
                case[_k_sec(k)][k] = v
            cases.append(case)

    return [dict(case) for case in cases]


if __name__ == "__main__":
    cases = config_cases(
        file_vars="../input/point_file_20220701.sfcf000.txt",
        infmt_opt=1,
        ntime=1,
        nlat=1,
        nlon=1,
        z0ghc=[0.001, 0.01],
        lambdars=[1.0, 1.25],
        product=True,
    )
    ds = run_config_sens(cases)
