DICT_CIS_VULNERABILITY_FLOOD = {
    
    "roads": {
                "motorway" : [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "motorway_link" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "trunk" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "trunk_link" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "primary" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "primary_link":  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "secondary" : ['F7.9','F7.8' ],
                "secondary_link" :  ['F7.9','F7.8' ],
                "tertiary" :  ['F7.9','F7.8' ],
                "tertiary_link" :  ['F7.9','F7.8' ],
                "residential" :  ['F7.9','F7.8' ],
                "road" :  ['F7.9','F7.8' ],
                "unclassified" :  ['F7.9','F7.8' ],
                "track" :  ['F7.9','F7.8' ],
                "service" :  ['F7.9','F7.8' ],},
    

    "main_roads" : {
                "motorway" : [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "motorway_link" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "trunk" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "trunk_link" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "primary" :  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "primary_link":  [ 'F7.5', 'F7.6', 'F7.7','F7.4',],
                "secondary" : ['F7.9','F7.8' ],
                "secondary_link" :  ['F7.9','F7.8' ],
                "tertiary" :  ['F7.9','F7.8' ],
                "tertiary_link" :  ['F7.9','F7.8' ],

        },

    "rail": {
        'rail' : ['F8.1', 'F8.2', 'F8.3', 'F8.4', 'F8.5', 'F8.6', 'F8.7'],
    },

    "air": { 
        "aerodrome" : ['F9.1','F9.2', 'F9.3'],
        "terminal" : ['F9.1','F9.2', 'F9.3'],
        "runway" : ['F7.4', 'F7.5', 'F7.6', 'F7.7'],},

    "telecom": { # curves beschikbaar voor overstromingen
        "mast" : ['F10.1'],
        "communications_tower" : ['F6.1', 'F6.2'],
        },

    "water_supply": { 
                "water_works" : ['F14.1','F14.2','F14.3','F14.4','F14.5','F14.6','F14.7','F14.8','F14.9','F14.10'], # curves beschikbaar voor overstromingen
                "water_well": ['F15.1'], # curves beschikbaar voor overstromingen
                "water_tower": ['F13.4'], 
                "reservoir_covered": ['F13.1', 'F13.2', 'F13.3', 'F13.4', 'F13.5'],
                "storage_tank": ['F13.1', 'F13.2', 'F13.3', 'F13.5'], # curves beschikbaar voor overstromingen
        },

    "waste_solid": {
            "waste_transfer_station" : ['F18.1'] 
            },

    "education": { # curves beschikbaar voor overstromingen
            "school" : ['F21.6', 'F21.7', 'F21.8', 'F21.10', 'F21.11', 'F21.13'], 
            "kindergarten" : ['F21.6', 'F21.7', 'F21.8', 'F21.10', 'F21.11', 'F21.13'], 
            "college": ['F21.6', 'F21.7', 'F21.8', 'F21.10', 'F21.11', 'F21.13'], 
            "university": ['F21.6', 'F21.7', 'F21.8', 'F21.10', 'F21.11', 'F21.13'], 
            "library": ['F21.6', 'F21.7', 'F21.8', 'F21.10', 'F21.11', 'F21.13'],


        },

    "healthcare": { # curves beschikbaar voor overstromingen
                "hospital" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'], 
                "clinic" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'], 
                "doctors" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "pharmacy" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "dentist" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "physiotherapist" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "alternative" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "laboratory" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "optometrist" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "rehabilitation" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "blood_donation" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
                "birthing_center" : ['F21.6', 'F21.8', 'F21.9', 'F21.12'],
        },

    "power": {
                "line" : ['F6.1', 'F6.2'], 
                "cable" : ['F5.1'], 
                "minor_line" : ['F6.1', 'F6.2'], 
                "plant" : ['F1.1', 'F1.2', 'F1.3', 'F1.4', 'F1.5', 'F1.6', 'F1.7'], 
                "generator" : ['F2.1', 'F2.2', 'F2.3'], 
                "substation" : ['F2.1', 'F2.2', 'F2.3'],
                "transformer" : ['F2.1', 'F2.2', 'F2.3'], 
                "pole" : ['F6.1', 'F6.2'] , 
                "portal" : ['F2.1', 'F2.2', 'F2.3'],
                "tower" : ['F6.1', 'F6.2'] ,
                "terminal" : ['F2.1', 'F2.2', 'F2.3'] , 
                "switch" : ['F2.1', 'F2.2', 'F2.3'],
                "catenary_mast" : ['F10.1'],
        },

    "gas": {
        "pipeline" : ['F16.1','F16.2','F16.3'], 
        "storage_tank" : ['F13.1', 'F13.2', 'F13.3', 'F13.5'], # use water storage tanks curves now
        "substation" : ['F2.1', 'F2.2', 'F2.3'],
        },

    "oil": {
        "substation" : ['F2.1', 'F2.2', 'F2.3'],
        "pipeline" : ['F16.1','F16.2','F16.3'],  
        "petroleum_well": ['F15.1'], # curves beschikbaar voor overstromingen
        "oil_refinery" : ['F1.4'] # uses the curve of thermal power plants
        },

    "waste_water": {
            "wastewater_plant" : ['F18.1', 'F18.2', 'F18.3', 'F18.4', 'F18.5', 'F18.6'],
            "waste_transfer_station": ['F18.1'],
        },

    "buildings": {
                "yes" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "house" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "residential" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "detached" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "hut" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "industrial" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "shed" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
                "apartments" : ['F21.1', 'F21.2', 'F21.3', 'F21.4', 'F21.5'],
        },
}