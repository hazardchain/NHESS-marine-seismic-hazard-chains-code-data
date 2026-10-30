"""Node layers and the 39 -> 15 class mapping, copied verbatim from
marine_hazard_chains_reproduction.Rmd so both analyses share one taxonomy."""

LAYER_FINE = {
    "SEQ": "Trigger", "EQ": "Trigger", "AS": "Trigger", "VE": "Trigger",
    "SLS": "Natural", "CLS": "Natural", "LS": "Natural", "RO": "Natural", "LB": "Natural", "SA": "Natural",
    "LI": "Natural", "GU": "Natural", "SUB": "Natural", "MB": "Natural", "TS": "Natural", "CFL": "Natural",
    "FF": "Natural", "MED": "Natural", "WP": "Natural", "WD": "Natural", "DI": "Natural", "FA": "Natural",
    "PHF": "Engineering", "SDF": "Engineering", "SPF": "Engineering", "SCF": "Engineering", "NFA": "Engineering",
    "FI": "Engineering", "EX": "Engineering", "CF": "Engineering", "NF": "Engineering",
    "MTD": "Consequence", "BI": "Consequence", "EC": "Consequence", "PP": "Consequence",
    "HD": "Consequence", "IR": "Consequence", "SU": "Consequence", "CP": "Consequence",
}

COARSE_MAP = {
    "SEQ": "EQ or AS", "AS": "EQ or AS", "EQ": "EQ or AS", "VE": "EQ or AS",
    "GU": "Mass Slide", "SLS": "Mass Slide", "CLS": "Mass Slide", "LI": "Mass Slide", "MB": "Mass Slide",
    "SUB": "Mass Slide", "LS": "Mass Slide", "RO": "Mass Slide", "LB": "Mass Slide", "SA": "Mass Slide",
    "PHF": "Critical Infra. Fail.", "SDF": "Critical Infra. Fail.", "CF": "Critical Infra. Fail.",
    "NFA": "Critical Infra. Fail.", "SPF": "Critical Infra. Fail.",
    "BI": "Business Interr.", "MTD": "Business Interr.", "WP": "Water Dil.", "WD": "Water Dil.",
    "NF": "Critical NET Fail.", "SCF": "Critical NET Fail.", "FA": "Flora an Fauna", "MED": "Flora an Fauna",
    "FI": "Fire", "EX": "Fire", "IR": "Healthcare Deg.", "HD": "Healthcare Deg.",
    "TS": "Flood", "CFL": "Flood", "FF": "Flood",
    "EC": "Econ. Crisis", "DI": "Disease", "PP": "Phy.-Psych. Trauma", "CP": "Culture Perils", "SU": "Social Unrest",
}

LAYER15 = {
    "EQ or AS": "Trigger",
    "Mass Slide": "Natural", "Water Dil.": "Natural", "Flora an Fauna": "Natural", "Flood": "Natural", "Disease": "Natural",
    "Critical Infra. Fail.": "Engineering", "Critical NET Fail.": "Engineering", "Fire": "Engineering",
    "Business Interr.": "Consequence", "Healthcare Deg.": "Consequence", "Econ. Crisis": "Consequence",
    "Phy.-Psych. Trauma": "Consequence", "Culture Perils": "Consequence", "Social Unrest": "Consequence",
}

ABBR15 = {
    "EQ or AS": "EQ", "Mass Slide": "MS", "Critical Infra. Fail.": "CIF", "Business Interr.": "BI",
    "Water Dil.": "WD", "Critical NET Fail.": "CNF", "Flora an Fauna": "FA", "Fire": "FI",
    "Healthcare Deg.": "HD", "Flood": "FL", "Econ. Crisis": "EC", "Disease": "DI",
    "Phy.-Psych. Trauma": "PP", "Culture Perils": "CP", "Social Unrest": "SU",
}

FULL_NAME = {
    "SEQ": "Submarine Earthquake", "EQ": "Coastal Earthquake", "AS": "Aftershock", "VE": "Volcanic Eruption",
    "SLS": "Submarine Landslide", "CLS": "Coastal Landslide", "LS": "Landslide", "RO": "Rockfall", "LB": "Lake Burst",
    "SA": "Snow Avalanche", "LI": "Liquefaction", "GU": "Ground Uplift", "SUB": "Subsidence", "MB": "Methane Blowout",
    "TS": "Tsunami", "CFL": "Coastal Flooding", "FF": "Fluvial Flood", "MED": "Marine Ecosystem Damage",
    "WP": "Water Pollution", "WD": "Water Dilemma", "DI": "Disease", "FA": "Flora and Fauna",
    "PHF": "Port / Harbour Failure", "SDF": "Seawall / Dike Failure", "SPF": "Subsea Pipeline Failure",
    "SCF": "Subsea Cable Failure", "NFA": "Nuclear Facility Accident", "FI": "Fire", "EX": "Explosion",
    "CF": "Critical Infrastructure Failure", "NF": "Critical Network Failure",
    "MTD": "Maritime Transport Disruption", "BI": "Business Interruption", "EC": "Economic Crisis",
    "PP": "Physical & Psychological Trauma", "HD": "Healthcare Degradation", "IR": "Insufficient Rescue",
    "SU": "Social Unrest", "CP": "Culture Perils",
}
