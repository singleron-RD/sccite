import json
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger("multiqc")
ASSAY = "sccite"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name=ASSAY,
            anchor=ASSAY,
        )
        log.info(f"Running module: {ASSAY}")

        stat_data = self.parse_json(ASSAY, "stats")
        if all(len(x) == 0 for x in [stat_data]):
            raise ModuleNoSamplesFound

        # Basic Stats Table
        self.general_stats_table(stat_data)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_json(self, assay, seg):
        data_dict = defaultdict(dict)
        n = 0
        for f in self.find_log_files(f"{assay}/{seg}"):
            log.debug(f"Found file: {f['fn']}")
            n += 1
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(f".{assay}")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {n} {assay} {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_{assay}_{seg}")
        return data_dict

    def general_stats_table(self, summary_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "Predefined pattern of barcode and UMI",
                "scale": "purple",
                "hidden": True,
            },
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Number of reads in the input file",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "Percent of reads with valid barcodes. Ideally > 80%",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Corrected Barcodes": {
                "title": "Corrected Barcodes",
                "description": "Percent of corrected barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Number of Cells": {
                "title": "N Cells",
                "description": "Number of cells from matched scRNA",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Fraction Tag Reads": {
                "title": "% Tag Reads",
                "description": "Fraction of valid reads that contain a recognized tag barcode. Ideally > 80%",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Fraction Tag Reads in Cell": {
                "title": "% Tag Reads in Cell",
                "description": "Fraction of tag reads that match the cell barcodes of scRNA. Ideally > 50%",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Median UMI per Cell": {
                "title": "Median UMI",
                "description": "Median number of tag UMI per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Aggregate Barcodes": {
                "title": "Aggregate Barcodes",
                "description": "Number of cell barcodes with protein aggregation. Ideally < 10",
                "scale": "blue",
                "format": "{:,.0f}",
            },
        }
        self.general_stats_addcols(summary_data, headers=headers)

    def barcode_rank_plot(self, umi_count_data):
        plot_data = {}
        colors = {}
        for sample in umi_count_data:
            for sub in umi_count_data[sample]:
                cur = umi_count_data[sample][sub]
                if not cur:
                    continue
                new = {}
                for k, v in cur.items():
                    new[int(k)] = v
                plot_data[sub] = new
                if "pure" in sub:
                    colors[sub] = "darkblue"
                elif "mix" in sub:
                    colors[sub] = "lightblue"
                elif "background" in sub:
                    colors[sub] = "lightgray"

        # Config for the plot
        pconfig = {
            "id": f"{ASSAY}_barcode_rank_plot",
            "title": f"{ASSAY}: Barcode Rank",
            "ylab": "Fragment Count",
            "xlab": "Barcode Rank",
            "yLog": True,
            "xLog": True,
            "colors": colors,
            "ymin": 0,
            "height": 750,
        }

        return linegraph.plot(plot_data, pconfig)

    def saturation_plot(self, saturation_data):
        # Config for the plot
        pconfig = {
            "id": "scrna_saturation_plot",
            "title": "scrna: Saturation",
            "ylab": "Saturation",
            "xlab": "Percent of Reads",
            "height": 750,
        }

        return linegraph.plot(saturation_data, pconfig)

    def median_gene_plot(self, median_gene_data):
        # Config for the plot
        pconfig = {
            "id": "scrna_median_gene_plot",
            "title": "scrna: Median Gene",
            "ylab": "Median Gene",
            "xlab": "Percent of Reads",
            "height": 750,
        }
        return linegraph.plot(median_gene_data, pconfig)
