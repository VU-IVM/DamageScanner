"""This file is part of OSM-flex.
Copyright (C) 2023 OSM-flex contributors listed in AUTHORS.
OSM-flex is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.
OSM-flex is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
-----
This file defines paths, URLs, and dictionaries used to:
- Manage OSM `.pbf` downloads via Geofabrik.
- Lookup region names from ISO3 country codes.
- Link infrastructure object types to predefined vulnerability curve sets.

Attributes:
    OSM_DIR (Path): Base directory for OSM data on the local machine.
    OSM_DATA_DIR (Path): Folder path where `.osm.pbf` files are stored.
    GEOFABRIK_URL (str): Base URL for Geofabrik region downloads.
    PLANET_URL (str): URL for downloading the full OpenStreetMap planet file.

Dictionaries:
    DICT_GEOFABRIK (dict): Maps ISO3 codes to Geofabrik region and subfolder.
    DICT_CIS_VULNERABILITY_FLOOD (dict): Maps infrastructure types to vulnerability curve IDs.
"""

from pathlib import Path

# =============================================================================
# FILE AND FOLDER PATHS
# =============================================================================
OSM_DIR = Path(Path.home(), "osm")
OSM_DATA_DIR = OSM_DIR.joinpath("osm_pbf")

# =============================================================================
# URLS
# =============================================================================

GEOFABRIK_URL = "https://download.geofabrik.de/"
PLANET_URL = "https://planet.openstreetmap.org/pbf/planet-latest.osm.pbf"


# =============================================================================
# DICTIONARIES
# =============================================================================

"""
dictionary for the generation of the correct download api-address at
geofabrik.de, relating ISO3-country codes to the region & written-out name.
Adapted from the GitHub repo osm_clipper (https://github.com/ElcoK/osm_clipper)
Used by download.get_country_geofabrik().

Note: A few small countries will be downloaded as a multi-country file, as
indicated in the comments.

Note: "special" ISO-3 codes - Canary Islands (IC), Asian part of Russia (RUS-A),
European part of Russia (RUS-E)
"""
DICT_GEOFABRIK = {
    "AFG": ("asia", "afghanistan"),
    "ALB": ("europe", "albania"),
    "DZA": ("africa", "algeria"),
    "AND": ("europe", "andorra"),
    "AGO": ("africa", "angola"),
    "BEN": ("africa", "benin"),
    "BWA": ("africa", "botswana"),
    "BFA": ("africa", "burkina-faso"),
    "BDI": ("africa", "burundi"),
    "CMR": ("africa", "cameroon"),
    "IC": ("africa", "canary-islands"),
    "CPV": ("africa", "cape-verde"),
    "CAF": ("africa", "central-african-republic"),
    "TCD": ("africa", "chad"),
    "COM": ("africa", "comores"),
    "COG": ("africa", "congo-brazzaville"),
    "COD": ("africa", "congo-democratic-republic"),
    "DJI": ("africa", "djibouti"),
    "EGY": ("africa", "egypt"),
    "GNQ": ("africa", "equatorial-guinea"),
    "ERI": ("africa", "eritrea"),
    "ETH": ("africa", "ethiopia"),
    "GAB": ("africa", "gabon"),
    "GMB": ("africa", "senegal-and-gambia"),  # TOGETHER WITH SENEGAL
    "GHA": ("africa", "ghana"),
    "GIN": ("africa", "guinea"),
    "GNB": ("africa", "guinea-bissau"),
    "CIV": ("africa", "ivory-coast"),
    "KEN": ("africa", "kenya"),
    "LSO": ("africa", "lesotho"),
    "LBR": ("africa", "liberia"),
    "LBY": ("africa", "libya"),
    "MDG": ("africa", "madagascar"),
    "MWI": ("africa", "malawi"),
    "MLI": ("africa", "mali"),
    "MRT": ("africa", "mauritania"),
    "MAR": ("africa", "morocco"),
    "ESH": ("africa", "morocco"),  # DOES INCLUDE THE WESTERN SAHARA
    "MOZ": ("africa", "mozambique"),
    "NAM": ("africa", "namibia"),
    "NER": ("africa", "niger"),
    "NGA": ("africa", "nigeria"),
    "RWA": ("africa", "rwanda"),
    "SHN": ("africa", "saint-helena-ascension-and-tristan-da-cunha"),
    "STP": ("africa", "sao-tome-and-principe"),
    "SEN": ("africa", "senegal-and-gambia"),  # TOGETHER WITH THE GAMBIA
    "SYC": ("africa", "seychelles"),
    "SLE": ("africa", "sierra-leone"),
    "SOM": ("africa", "somalia"),  # TOGETHER WITH SOMALILAND
    "SOL": ("africa", "somaliland"),  # TOGETHER WITH SOMALIA
    "ZAF": ("africa", "south-africa"),
    "SDN": ("africa", "sudan"),
    "SSD": ("africa", "south-sudan"),
    "SWZ": ("africa", "swaziland"),
    "TZA": ("africa", "tanzania"),
    "TGO": ("africa", "togo"),
    "TUN": ("africa", "tunisia"),
    "UGA": ("africa", "uganda"),
    "ZMB": ("africa", "zambia"),
    "ZWE": ("africa", "zimbabwe"),
    "MUS": ("africa", "mauritius"),
    "ARM": ("asia", "armenia"),
    "AZE": ("asia", "azerbaijan"),
    "BGD": ("asia", "bangladesh"),
    "BTN": ("asia", "bhutan"),
    "KHM": ("asia", "cambodia"),
    "CHN": ("asia", "china"),
    "SAU": (
        "asia",
        "gcc-states",
    ),  # Together with Kuwait, the United Arab Emirates, Qatar, Bahrain, and Oman
    "KWT": (
        "asia",
        "gcc-states",
    ),  # Together with Saudi Arabia, the United Arab Emirates, Qatar, Bahrain, and Oman
    "ARE": (
        "asia",
        "gcc-states",
    ),  # Together with Saudi Arabia, Kuwait, Qatar, Bahrain, and Oman
    "QAT": (
        "asia",
        "gcc-states",
    ),  # Together with Saudi Arabia, Kuwait, the United Arab Emirates, Bahrain, and Oman
    "OMN": (
        "asia",
        "gcc-states",
    ),  # Together with Saudi Arabia, Kuwait, the United Arab Emirates, Qatar and Oman
    "BHR": (
        "asia",
        "gcc-states",
    ),  # Together with Saudi Arabia, Kuwait, the United Arab Emirates, Qatar and Bahrain
    "IND": ("asia", "india"),
    "IDN": ("asia", "indonesia"),
    "IRN": ("asia", "iran"),
    "IRQ": ("asia", "iraq"),
    "ISR": ("asia", "israel-and-palestine"),  # TOGETHER WITH PALESTINE
    "PSE": ("asia", "israel-and-palestine"),  # TOGETHER WITH ISRAEL
    "JPN": ("asia", "japan"),
    "JOR": ("asia", "jordan"),
    "KAZ": ("asia", "kazakhstan"),
    "KGZ": ("asia", "kyrgyzstan"),
    "LAO": ("asia", "laos"),
    "LBN": ("asia", "lebanon"),
    "MYS": ("asia", "malaysia-singapore-brunei"),  # TOGETHER WITH SINGAPORE AND BRUNEI
    "SGP": ("asia", "malaysia-singapore-brunei"),  # TOGETHER WITH MALAYSIA AND BRUNEI
    "BRN": (
        "asia",
        "malaysia-singapore-brunei",
    ),  # TOGETHER WITH MALAYSIA AND SINGAPORE
    "MDV": ("asia", "maldives"),
    "MNG": ("asia", "mongolia"),
    "MMR": ("asia", "myanmar"),
    "NPL": ("asia", "nepal"),
    "PRK": ("asia", "north-korea"),
    "PAK": ("asia", "pakistan"),
    "PHL": ("asia", "philippines"),
    "RUS-A": ("asia", "russia"),  # Asian part of Russia
    "KOR": ("asia", "south-korea"),
    "LKA": ("asia", "sri-lanka"),
    "SYR": ("asia", "syria"),
    "TWN": ("asia", "taiwan"),
    "TJK": ("asia", "tajikistan"),
    "THA": ("asia", "thailand"),
    "TLS": ("asia", "east-timor"),
    "TKM": ("asia", "turkmenistan"),
    "UZB": ("asia", "uzbekistan"),
    "VNM": ("asia", "vietnam"),
    "YEM": ("asia", "yemen"),
    "BHS": ("central-america", "bahamas"),
    "BLZ": ("central-america", "belize"),
    "CUB": ("central-america", "cuba"),
    "GTM": ("central-america", "guatemala"),
    "HTI": ("central-america", "haiti-and-domrep"),  # TOGETHER WITH DOMINICAN REPUBLIC
    "DOM": ("central-america", "haiti-and-domrep"),  # TOGETHER WITH HAITI
    "JAM": ("central-america", "jamaica"),
    "HND": ("central-america", "honduras"),
    "NIC": ("central-america", "nicaragua"),
    "SLV": ("central-america", "el-salvador"),
    "CRI": ("central-america", "costa-rica"),
    "PAN": ("central-america", "panama"),
    "AUT": ("europe", "austria"),
    "BLR": ("europe", "belarus"),
    "BEL": ("europe", "belgium"),
    "BIH": ("europe", "bosnia-herzegovina"),
    "BGR": ("europe", "bulgaria"),
    "HRV": ("europe", "croatia"),
    "CYP": ("europe", "cyprus"),
    "XAD": ("europe", "cyprus"),  # DOES INCLUDE AKROTIRI AND DHEKELIA
    "ZNC": ("europe", "cyprus"),  # DOES INCLUDE NORTHERN CYPRUS
    "CZE": ("europe", "czech-republic"),
    "DNK": ("europe", "denmark"),
    "EST": ("europe", "estonia"),
    "FRO": ("europe", "faroe-islands"),
    "FIN": ("europe", "finland"),
    "ALA": ("europe", "finland"),  # TOGETHER WITH ALAND
    "FRA": ("europe", "france"),
    "GEO": ("europe", "georgia"),
    "DEU": ("europe", "germany"),
    "GBR": ("europe", "great-britain"),  # DOES NOT INCLUDE NORTHERN ISLAND
    "GRC": ("europe", "greece"),
    "HUN": ("europe", "hungary"),
    "ISL": ("europe", "iceland"),
    "IRL": ("europe", "ireland-and-northern-ireland"),
    "IMN": ("europe", "isle-of-man"),
    "ITA": ("europe", "italy"),
    "LVA": ("europe", "latvia"),
    "LIE": ("europe", "liechtenstein"),
    "LTU": ("europe", "lithuania"),
    "LUX": ("europe", "luxembourg"),
    "MKD": ("europe", "macedonia"),
    "MLT": ("europe", "malta"),
    "MDA": ("europe", "moldova"),
    "MCO": ("europe", "monaco"),
    "MNE": ("europe", "montenegro"),
    "NLD": ("europe", "netherlands"),
    "NOR": ("europe", "norway"),
    "SJM": ("europe", "norway"),  # DOES INCLUDE SVALBARD AND JAN MAYEN
    "POL": ("europe", "poland"),
    "PRT": ("europe", "portugal"),
    "ROU": ("europe", "romania"),
    "RUS-E": ("europe", "russia"),  # European part of Russia
    "SRB": ("europe", "serbia"),
    "SVK": ("europe", "slovakia"),
    "SVN": ("europe", "slovenia"),
    "ESP": ("europe", "spain"),
    "SWE": ("europe", "sweden"),
    "CHE": ("europe", "switzerland"),
    "TUR": ("europe", "turkey"),
    "UKR": ("europe", "ukraine"),
    "XKO": ("europe", "kosovo"),
    "SMR": ("europe", "italy"),  # DOES INCLUDE SAN MARINO
    "VAT": ("europe", "italy"),  # DOES INCLUDE VATICAN CITY
    "GGY": ("europe", "guernsey-jersey"),  # TOGETHER WITH JERSEY
    "JEY": ("europe", "guernsey-jersey"),  # TOGETHER WITH GUERNSEY
    "CAN": ("north-america", "canada"),
    "SPM": ("north-america", "canada"),  # DOES INCLUDE SAINT PIERRE AND MIQUELON
    "GRL": ("north-america", "greenland"),
    "MEX": ("north-america", "mexico"),
    "USA": ("north-america", "us"),
    "AUS": ("australia-oceania", "australia"),
    "CXR": ("australia-oceania", "australia"),  # DOES INCLUDE CHRISTMAS ISLAND
    "CCK": ("australia-oceania", "australia"),  # DOES INCLUDE COCOS ISLAND
    "HMD": (
        "australia-oceania",
        "australia",
    ),  # DOES INCLUDE HEARD ISLAND AND MCDONALD ISLAND
    "NFK": ("australia-oceania", "australia"),  # DOES INCLUDE HEARD NORFOLK ISLAND
    "COK": ("australia-oceania", "cook-islands"),
    "FJI": ("australia-oceania", "fiji"),
    "KIR": ("australia-oceania", "kiribati"),
    "MHL": ("australia-oceania", "marshall-islands"),
    "FSM": ("australia-oceania", "micronesia"),
    "NRU": ("australia-oceania", "nauru"),
    "NCL": ("australia-oceania", "new-caledonia"),
    "NZL": ("australia-oceania", "new-zealand"),
    "NIU": ("australia-oceania", "niue"),
    "PLW": ("australia-oceania", "palau"),
    "PNG": ("australia-oceania", "papua-new-guinea"),
    "PCN": ("australia-oceania", "pitcairn-islands"),
    "WSM": ("australia-oceania", "samoa"),
    "SLB": ("australia-oceania", "solomon-islands"),
    "TON": ("australia-oceania", "tonga"),
    "TUV": ("australia-oceania", "tuvalu"),
    "VUT": ("australia-oceania", "vanuatu"),
    "PYF": ("australia-oceania", "polynesie-francaise"),
    "GUM": ("australia-oceania", "american-oceania"),  # DOES INCLUDE GUAM
    "MNP": (
        "australia-oceania",
        "american-oceania",
    ),  # DOES INCLUDE NORTHERN MARIANA ISLANDS
    "ARG": ("south-america", "argentina"),
    "BOL": ("south-america", "bolivia"),
    "BRA": ("south-america", "brazil"),
    "CHL": ("south-america", "chile"),
    "COL": ("south-america", "colombia"),
    "ECU": ("south-america", "ecuador"),
    "GUY": ("south-america", "guyana"),
    "PRY": ("south-america", "paraguay"),
    "PER": ("south-america", "peru"),
    "SUR": ("south-america", "suriname"),
    "URY": ("south-america", "uruguay"),
    "VEN": ("south-america", "venezuela"),
}


DICT_CIS_VULNERABILITY_FLOOD = {
    "roads": {
        "motorway": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "motorway_link": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "trunk": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "trunk_link": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "primary": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "primary_link": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "secondary": ["F7.9", "F7.8"],
        "secondary_link": ["F7.9", "F7.8"],
        "tertiary": ["F7.9", "F7.8"],
        "tertiary_link": ["F7.9", "F7.8"],
        "residential": ["F7.9", "F7.8"],
        "road": ["F7.9", "F7.8"],
        "unclassified": ["F7.9", "F7.8"],
        "track": ["F7.9", "F7.8"],
        "service": ["F7.9", "F7.8"],
    },
    "main_roads": {
        "motorway": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "motorway_link": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "trunk": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "trunk_link": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "primary": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "primary_link": [
            "F7.5",
            "F7.6",
            "F7.7",
            "F7.4",
        ],
        "secondary": ["F7.9", "F7.8"],
        "secondary_link": ["F7.9", "F7.8"],
        "tertiary": ["F7.9", "F7.8"],
        "tertiary_link": ["F7.9", "F7.8"],
    },
    "rail": {
        "rail": ["F8.1", "F8.2", "F8.3", "F8.4", "F8.5", "F8.6", "F8.7"],
    },
    "air": {
        "aerodrome": ["F9.1", "F9.2", "F9.3"],
        "terminal": ["F9.1", "F9.2", "F9.3"],
        "runway": ["F7.4", "F7.5", "F7.6", "F7.7"],
    },
    "telecom": {  # curves beschikbaar voor overstromingen
        "mast": ["F10.1"],
        "tower": ["F6.1", "F6.2"],
    },
    "water_supply": {
        "water_works": [
            "F14.1",
            "F14.2",
            "F14.3",
            "F14.4",
            "F14.5",
            "F14.6",
            "F14.7",
            "F14.8",
            "F14.9",
            "F14.10",
        ],  # curves beschikbaar voor overstromingen
        "water_well": ["F15.1"],  # curves beschikbaar voor overstromingen
        "water_tower": ["F13.4"],
        "reservoir_covered": ["F13.1", "F13.2", "F13.3", "F13.4", "F13.5"],
        "storage_tank": [
            "F13.1",
            "F13.2",
            "F13.3",
            "F13.5",
        ],  # curves beschikbaar voor overstromingen
    },
    "waste_solid": {"waste_transfer_station": ["F18.1"]},
    "education": {  # curves beschikbaar voor overstromingen
        "school": ["F21.6", "F21.7", "F21.8", "F21.10", "F21.11", "F21.13"],
        "kindergarten": ["F21.6", "F21.7", "F21.8", "F21.10", "F21.11", "F21.13"],
        "college": ["F21.6", "F21.7", "F21.8", "F21.10", "F21.11", "F21.13"],
        "university": ["F21.6", "F21.7", "F21.8", "F21.10", "F21.11", "F21.13"],
        "library": ["F21.6", "F21.7", "F21.8", "F21.10", "F21.11", "F21.13"],
    },
    "healthcare": {  # curves beschikbaar voor overstromingen
        "hospital": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "clinic": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "doctors": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "pharmacy": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "dentist": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "physiotherapist": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "alternative": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "laboratory": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "optometrist": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "rehabilitation": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "blood_donation": ["F21.6", "F21.8", "F21.9", "F21.12"],
        "birthing_center": ["F21.6", "F21.8", "F21.9", "F21.12"],
    },
    "power": {
        "line": ["F6.1", "F6.2"],
        "cable": ["F5.1"],
        "minor_line": ["F6.1", "F6.2"],
        "plant": ["F1.1", "F1.2", "F1.3", "F1.4", "F1.5", "F1.6", "F1.7"],
        "generator": ["F2.1", "F2.2", "F2.3"],
        "substation": ["F2.1", "F2.2", "F2.3"],
        "transformer": ["F2.1", "F2.2", "F2.3"],
        "pole": ["F6.1", "F6.2"],
        "portal": ["F2.1", "F2.2", "F2.3"],
        "tower": ["F6.1", "F6.2"],
        "terminal": ["F2.1", "F2.2", "F2.3"],
        "switch": ["F2.1", "F2.2", "F2.3"],
        "catenary_mast": ["F10.1"],
    },
    "gas": {
        "pipeline": ["F16.1", "F16.2", "F16.3"],
        "storage_tank": [
            "F13.1",
            "F13.2",
            "F13.3",
            "F13.5",
        ],  # use water storage tanks curves now
        "substation": ["F2.1", "F2.2", "F2.3"],
    },
    "oil": {
        "substation": ["F2.1", "F2.2", "F2.3"],
        "pipeline": ["F16.1", "F16.2", "F16.3"],
        "petroleum_well": ["F15.1"],  # curves beschikbaar voor overstromingen
        "oil_refinery": ["F1.4"],  # uses the curve of thermal power plants
    },
    "waste_water": {
        "wastewater_plant": ["F18.1", "F18.2", "F18.3", "F18.4", "F18.5", "F18.6"],
        "waste_transfer_station": ["F18.1"],
    },
    "buildings": {
        "yes": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "house": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "residential": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "detached": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "hut": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "industrial": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "shed": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
        "apartments": ["F21.1", "F21.2", "F21.3", "F21.4", "F21.5"],
    },
}
