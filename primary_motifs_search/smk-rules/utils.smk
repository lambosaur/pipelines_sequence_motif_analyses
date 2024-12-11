from dotmap import DotMap
import os
from glom import glom


class DefaultDotMap(DotMap):
    """Default DotMap dict-with-dot-notation-access that returns None if key is missing."""

    def __getitem__(self, key):
        return super().__getitem__(key) if key in self else None


def prefix_with_basedir(fp: os.PathLike) -> str:
    """If provided a relative filepath, transform to an absolute filepath from its relation to the workflow.basedir."""
    if not os.path.isabs(fp):
        return os.path.normpath(Path(workflow.basedir) / fp)
    return str(fp)


def get_fasta_filepath(wildcards: Wildcards) -> str:
    fp = config["filepaths"][
        CONFIG_EXPERIMENTS[wildcards.EXPERIMENT_ID]["sequences_filepath_id"]
    ]
    return prefix_with_basedir(fp)


def get_output_dir() -> str:
    """Retrieve `config.filepaths.output_dir` or `config.filepaths.{PIPELINE_NAME}.output_dir` if defined, defaults to "results/{PIPELINE_NAME}".

    NOTE: the retrieved output dir path is prefixed with the workflow.basedir, if not an absolute path.
    """

    if "output_dir" in config["filepaths"]:
        # Generic output dir for all sub-pipelines in `primary_motifs_search`
        output_dir = Path(config["filepaths"]["output_dir"]) / PIPELINE_NAME

    elif (
        glom(
            target=config["filepaths"], spec=f"{PIPELINE_NAME}.output_dir", default=None
        )
        is not None
    ):
        # Specific output dir for this sub-pipeline
        output_dir = Path(config["filepaths"][PIPELINE_NAME]["output_dir"])

    else:
        # Default output dir
        output_dir = Path("results") / PIPELINE_NAME

    return prefix_with_basedir(output_dir)


def make_output_filepath(output_filepath: os.PathLike) -> str:
    """Return a normalized path of the output_filepath prefixed by retrieved output dir, if output_filepath is relative."""
    if os.path.isabs(output_filepath):
        return str(output_filepath)
    else:
        output_dir: Path = Path(get_output_dir())
        return str(os.path.normpath(output_dir / output_filepath))


def get_mandatory_filepath_from_experiments_config(
    wildcards: Wildcards, keychain: str
) -> str:
    """Retrieve a dotted key-chain from the `CONFIG_EXPERIMENTS.experiments` config ; raise KeyError if not found.

    NOTE: the first part of the keychain is expected to be an experiment_id.

    Raises KeyError if the keychain is not found in the experiments config.
    """
    # This will raise an error if not defined / found.
    try:
        FILEPATH_ID = glom(
            target=CONFIG_EXPERIMENTS.experiments,
            spec=keychain,
        )
    except PathAccessError:
        raise KeyError(
            f"Could not find key {keychain} in config CONFIG_EXPERIMENTS.experiments."
        )

    fp = config["filepaths"][FILEPATH_ID]
    return prefix_with_basedir(fp)


def get_optional_filepath_from_experiments_config(
    wildcards: Wildcards, keychain: str
) -> Union[str, List[Never]]:
    """Return filepath if keychain exists in experiments config, otherwise empty list (mandatory for Snakemake input)."""
    FILEPATH_ID = glom(
        target=CONFIG_EXPERIMENTS.experiments,
        spec=keychain,
        default=None,
    )
    if FILEPATH_ID is None:
        return []
    else:
        fp = glom(target=config["filepaths"], spec=FILEPATH_ID, default=None)
        if fp is None:
            return []
        return prefix_with_basedir(fp)
