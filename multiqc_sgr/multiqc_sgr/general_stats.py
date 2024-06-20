import json
import logging

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
log = logging.getLogger(__name__)


def get_frac(x):
    return round(float(x) * 100, 2)


def get_int(x):
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="scsnp",
            anchor="scsnp",
            info="Single cell amplicon variant calling.",
        )

        # Find and load any STAR reports
        count_data = self.parse_json("count")
        if all(len(x) == 0 for x in [count_data]):
            raise ModuleNoSamplesFound

        # assgin plot
        self.add_section(name="Variant count", anchor="scsnp_count", plot=self.count_bar(count_data))

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_json(self, seg):
        data_dict = dict()
        for f in self.find_log_files(f"scsnp/{seg}"):
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                s_name = f["s_name"].removesuffix(f".{seg}")
                if s_name in data_dict:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name] = parsed_data

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {len(data_dict)} scsnp {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_scsnp_{seg}")
        return data_dict


def general_stats_table(self, summary_data):
    protocol = {
        "protocol": {
            "title": "Protocol",
            "description": "Barcode and UMI protocol",
        },
        "Q30 of Barcodes": {
            "title": "Q30 of Barcodes",
            "description": "Q30 of barcodes",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "green",
            "modify": get_frac,
            "format": "{:,.2f}",
        },
    }
    starsolo = {
        "Valid Reads": {
            "title": "% Valid Reads",
            "description": "Fraction of reads with valid barcodes and UMI",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "green",
            "modify": get_frac,
            "format": "{:,.2f}",
        },
        "Estimated Number of Cells": {
            "title": "N Cells",
            "description": "Estimated number of cells",
            "format": get_int,
        },
        "Fraction Reads in Cells": {
            "title": "% Reads in Cells",
            "description": "Fraction of unique reads in cells",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "green",
            "modify": get_frac,
            "format": "{:,.2f}",
        },
        "Median Genes per Cell": {
            "title": "Median Genes",
            "description": "Median genes per cell",
            "format": get_int,
        },
        "Mean Used Reads per Cell": {
            "title": "Mean Reads",
            "description": "Mean used reads per cell",
            "format": get_int,
        },
        "Mean UMI per Cell": {
            "title": "Mean UMI",
            "description": "Mean UMI per Cell",
            "format": get_int,
        },
        "Saturation": {
            "title": "Saturation",
            "description": "Sequencing Saturation",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "green",
            "modify": get_frac,
            "format": "{:,.2f}",
        },
    }
    headers = {}
    for h in [protocol, starsolo]:
        headers.update(h)

    self.general_stats_addcols(summary_data, headers=headers)
