def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "sccite/stats": {"fn": "*sccite.*stats.json"},
    }
    config.update_dict(config.sp, sgr_search_patterns)
