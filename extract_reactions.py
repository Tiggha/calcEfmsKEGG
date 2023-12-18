import xml.etree.ElementTree as ET
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model

import os

from cobra.util.array import *

mod_name = ""
model = Model(mod_name + "_model")


def find_unique_entries(substrates, products):
    unique_substrates = set(substrates) - set(products)
    unique_products = set(products) - set(substrates)

    return list(unique_substrates), list(unique_products)


# Define a function to parse KGML and extract reactions, bounds, and stoichiometry
def parse_kgml(kgml_file, mod_name):
    reactions = {}
    reaction_bounds = {}
    stoichiometry = {}
    coefficient = 1
    lower_bound = 0
    upper_bound = 1000
    metabolite_dict = {}
    substrates = []
    products = []

    tree = ET.parse(kgml_file)
    root = tree.getroot()

    for entry in root.findall(".//reaction"):
        reaction_id = entry.get("id")
        reactions[reaction_id] = entry.get("name")
        reaction_name = entry.get("name")
        reaction_bounds[reaction_id] = entry.get("type")
        stoichiometry[reaction_id] = {}
        # print(entry)
        # print(reaction_id)
        # print(reaction_name.split(":", 1)[1])

        if entry.get("type") == "reversible":
            lower_bound = -1000
            upper_bound = 1000
        elif entry.get("type") == "irreversible":
            lower_bound = 0
            upper_bound = 1000
        rea_name_split = reaction_name.split(":", 1)[1]
        reaction = Reaction(rea_name_split.replace(" ", ""))
        reaction.name = reaction_name
        reaction.lower_bound = lower_bound
        reaction.upper_bound = upper_bound

        for substrate in entry.findall("substrate"):
            substrates.append(substrate.get("name"))
            metabolite_id = substrate.get("name")
            if " " in metabolite_id:
                metabolite_id = metabolite_id.split()[0]

            if metabolite_id not in metabolite_dict:
                metabolite = Metabolite(metabolite_id, compartment="default")
                metabolite_dict[metabolite_id] = metabolite
            else:
                metabolite = metabolite_dict[metabolite_id]
            model.add_metabolites(metabolite)
            reaction.add_metabolites({metabolite: -coefficient})

            # coefficient = float(substrate.get("stoichiometry"))
            stoichiometry[reaction_id][
                metabolite_id
            ] = -coefficient  # Substrates have negative coefficients

        for product in entry.findall("product"):
            products.append(product.get("name"))
            metabolite_id = product.get("name")
            if " " in metabolite_id:
                metabolite_id = metabolite_id.split()[0]

            if metabolite_id not in metabolite_dict:
                metabolite = Metabolite(metabolite_id, compartment="default")
                metabolite_dict[metabolite_id] = metabolite
            else:
                metabolite = metabolite_dict[metabolite_id]
            model.add_metabolites(metabolite)
            reaction.add_metabolites({metabolite: coefficient})

            # coefficient = float(product.get("stoichiometry"))
            stoichiometry[reaction_id][
                metabolite_id
            ] = coefficient  # Products have positive coefficients
        # print(reaction.reaction)

        model.add_reaction(reaction)
    unique_substrates, unique_products = find_unique_entries(substrates, products)

    unique_metabolites = unique_substrates + unique_products
    for i, met in enumerate(unique_metabolites):
        reaction = Reaction("R_EX_{}".format(i + 1))
        reaction.name = met
        reaction.lower_bound = -1000
        reaction.upper_bound = 1000

        meta = ""
        for x in model.metabolites:
            if x.id == met:
                meta = x
        try:
            reaction.add_metabolites({meta: -coefficient})
        except ValueError:
            print("error meta: ", meta)

        model.add_reaction(reaction)
    print("Unique Substrates:", unique_substrates)
    print("Unique Products:", unique_products)

    mypath = "extracted_KGML"
    with open(
        mypath + "/" + "efms/exchange/" + "exchange_" + mod_name + ".unique", "w"
    ) as f:
        f.write("exchange reactions: " + str(len(unique_metabolites)))

    return reactions, reaction_bounds, stoichiometry


np.set_printoptions(
    edgeitems=30,
    linewidth=100000,
)

# path of kgml files
mypath = "KGML"

# path of sbml files to be extracted to
extractPath = "extracted_KGML"
kgmlfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

for kgml in kgmlfiles:
    print("extracting: ", kgml)
    kgml_name = kgml.split(".")
    mod_name = kgml_name[0]
    model = Model(mod_name + "_model")

    reactions, reaction_bounds, stoichiometry = parse_kgml(
        mypath + "/" + kgml, mod_name
    )
    write_sbml_model(model, extractPath + "/" + mod_name + ".sbml.xml")

# Print the extracted information
""" for reaction_id, reaction_name in reactions.items():
    print(f"Reaction ID: {reaction_id}")
    print(f"Reaction Name: {reaction_name}")
    print(f"Bounds: {reaction_bounds[reaction_id]}")
    print("Stoichiometry:")
    for compound, coefficient in stoichiometry[reaction_id].items():
        print(f"{compound}: {coefficient}")
    print("\n")
print(model.reactions)
print(model.metabolites) """
