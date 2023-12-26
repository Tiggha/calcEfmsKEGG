import efmtool
import cobra
import subprocess
import distutils.spawn
import numpy
import requests
import os
import json
import libsbml


from pathlib import Path
from time import time
from guppy import hpy
from efmlrs import pre, post
from cobra.util.array import *
from shutil import copyfile
from argparse import ArgumentParser

from cobra.io import write_sbml_model
from cobra import Model, Reaction, Metabolite

KGML = "models/"
RESULT = "result/"
IGNORE_COMPARTIMENT = [""]
BOUNDS = False
LRS = "lrslib-071b"
KEGGtranslator = "KEGGtranslator_v2.5.jar"

model_name = ""
model_mfile = ""
model_rfile = ""
model_rvfile = ""
model_sfile = ""
efms_efmtool = ""
efms_mplrs = ""
info = ""
output_mplrs = ""
ine_cmp = ""


class efms:
    def __init__(self) -> None:
        self._reactions = []
        self._metabolites = []
        self._stoichiometry = None
        self._reversibilities = None
        self._efms = None

    def get_reactions(self):
        return len(self._reactions)

    def get_metabolites(self):
        return len(self._metabolites)

    def get_stoichiometry(self):
        return len(self._stoichiometry)

    def get_reversibilities(self):
        return np.count_nonzero(self._reversibilities)

    def get_irreversibilities(self):
        return np.count_nonzero(self._reversibilities == 0)

    def get_efms(self):
        return self._efms

    def set_reversibilityTrue():
        lower_bound = -1000
        upper_bound = 1000
        return lower_bound, upper_bound

    def set_reversibilityFalse():
        lower_bound = 0
        upper_bound = 1000
        return lower_bound, upper_bound

    def set_fileNames(self, isSBML, model_namee):
        global model_name
        global model_mfile
        global model_rfile
        global model_rvfile
        global model_sfile
        global efms_efmtool
        global efms_mplrs
        global info
        global output_mplrs
        global ine_cmp

        model_name = model_namee

        if isSBML:
            model_mfile = model_name + ".sbml.mfile"
            model_rfile = model_name + ".sbml.rfile"
            model_rvfile = model_name + ".sbml.rvfile"
            model_sfile = model_name + ".sbml.sfile"
            efms_efmtool = model_name + ".sbml_cmp_efmtool.efms"
            efms_mplrs = model_name + ".sbml_cmp_mplrs.efms"
            info = model_name + ".sbml.info"
            output_mplrs = model_name + ".sbml_decomp_mplrs.efms"
            ine_cmp = KGML + model_name + "/" + model_name + ".sbml_cmp.ine"

        else:
            model_mfile = model_name + ".mfile"
            model_rfile = model_name + ".rfile"
            model_rvfile = model_name + ".rvfile"
            model_sfile = model_name + ".sfile"
            efms_efmtool = model_name + "_cmp_efmtool.efms"
            efms_mplrs = model_name + "_cmp_mplrs.efms"
            info = model_name + ".info"
            output_mplrs = model_name + "_decomp_mplrs.efms"
            ine_cmp = KGML + model_name + "/" + model_name + "_cmp.ine"

    def _safety_checks(self):
        """
        Check if `mpirun` is installed (executable found?).
        Raises
        ------
        OSError
            Raises error if executable `mpirun` cannot be found.
        """
        mpirun = distutils.spawn.find_executable("mpirun")

        if mpirun is None:
            raise OSError("mpirun cannot be located. Is it installed?")

    def _checkForReactions(self, model):
        if not os.path.isfile(model):
            "File/Model does not exist."
            exit()

        with open(model) as f:
            if "reaction" in f.read():
                print("Model has reactions.")
            else:
                print("Model has no reactions.")
                exit()

    def _readReversiblitiesFromSbml(self, sbml_file):
        reactions_reversibilities = {}
        reader = libsbml.SBMLReader()
        sbml = reader.readSBML(sbml_file)
        model = sbml.getModel()
        for r in range(model.getNumReactions()):
            rxn = model.getReaction(r)
            rev = rxn.reversible
            reaction_id = rxn.getId()
            reactions_reversibilities[reaction_id] = rev
        return reactions_reversibilities

    def _get_sbml(self, sbml_file: str) -> None:
        """
        get model, stoichiometry, reversibilities,
        reactions and metabolites from file
        """
        sbml_model = cobra.io.read_sbml_model(sbml_file)
        self._stoichiometry = cobra.util.array.create_stoichiometric_matrix(
            sbml_model
        )  # noqa

        reactions_reversibilities = self._readReversiblitiesFromSbml(sbml_file)
        if hasattr(sbml_model.reactions[0], "reversible"):
            # KEGG models
            self._reversibilities = np.array(
                [rea.reversible for rea in sbml_model.reactions]
            ).astype(int)

        else:
            # non-KEGG models

            self._reversibilities = np.array(
                [rea.reversibility for rea in sbml_model.reactions]
            ).astype(int)

        for x in sbml_model.reactions:
            self._reactions.append(x.id)
        for x in sbml_model.metabolites:
            self._metabolites.append(x.id)

        return (
            self._stoichiometry,
            self._reversibilities,
            self._reactions,
            self._metabolites,
        )

    def _from_KEGG(self, model):
        self._getKGML(model)
        self._runKeggTranslatorToSbml()
        return self._get_sbml(KGML + model + ".sbml.xml")

    def _calc_efms_mplrs(self, path_to_lrs, infile, outfile):
        """
        calculate efms with mplrs
        """
        subprocess.check_output(["mpirun", path_to_lrs + "/mplrs", infile, outfile])

    def supp(self, vector):
        return list(np.nonzero(np.round(vector, 10))[0])

    def _calc_efms_efmtool(self, model) -> None:
        """
        calculate efms with efmtool
        """
        self._get_sbml(model)

        np.set_printoptions(threshold=np.inf)

        efms = efmtool.calculate_efms(
            self._stoichiometry,
            self._reversibilities,
            self._reactions,
            self._metabolites,
            options=None,
            jvm_options=None,
        )
        efm_supps = [tuple(self.supp(efm)) for efm in efms.T]

        unique_efm_supps = set(efm_supps)

        self._efms = len(unique_efm_supps)
        numpy.set_printoptions(suppress=True)
        numpy.savetxt(
            KGML + model_name + "/" + efms_efmtool, efms, fmt="%1.1f", delimiter="\t"
        )

        return efms

    def _create_custom_model(
        self, name_model, stoichiometry, reversibilities, reactions, metabolites
    ):
        """
        stoichiometry: np.array of numbers
        reversibilities: array of strings
        reactions: array of strings
        metabolites: array of strings
        """
        # Define the stoichiometric matrix
        S = stoichiometry  # np.array([[1, -1, 0, -1, 0], [0, 1, -1, 0, 0], [0, 0, 1, 1, -1]])

        # Define the reversibility
        rev = reversibilities  # [0, 1, 1, 0, 1]

        # Create a new cobrapy model
        model = Model(name_model)

        # Add the metabolites to the model
        metabolites = 0
        for i, met in enumerate(S):
            metabolite = Metabolite("M_{}".format(i + 1), compartment="c")
            model.add_metabolites(metabolite)
            metabolites = i + 1

        # Add the reactions to the model
        for i, col in enumerate(S.T):
            reaction = Reaction("R_{}".format(i + 1))
            reversibility = bool(rev[i])
            if reversibility == True:
                reaction.lower_bound = self.set_reversibilityTrue()[0]
                reaction.upper_bound = self.set_reversibilityTrue()[1]
            else:
                reaction.lower_bound = self.set_reversibilityFalse()[0]
                reaction.upper_bound = self.set_reversibilityFalse()[1]
            reaction.add_metabolites(
                {model.metabolites[j]: float(S[j, i]) for j in range(metabolites)}
            )  # noqa
            model.add_reaction(reaction)

        # Save the model to an XML file
        write_sbml_model(model, name_model + ".xml")

    def _get_model_name(self, model):
        filename = os.path.basename(model)
        filename2 = os.path.splitext(filename)
        model_name = os.path.splitext(filename2[0])[0]
        return model_name

    def _calc_efms_lrs(self):
        self._preprocessing()
        self._calc_efms_mplrs(LRS, ine_cmp, KGML + model_name + "/" + efms_mplrs)
        self._postprocessing_mplrs()
        efms = self._get_efms_from_file()
        return efms

    def _create_dir_for_model(self):
        Path(KGML + model_name).mkdir(parents=True, exist_ok=True)
        # copy model
        if os.path.isfile(KGML + model_name + ".sbml.xml"):
            copyfile(
                KGML + model_name + ".sbml.xml",
                KGML + model_name + "/" + model_name + ".sbml.xml",
            )

        if os.path.isfile(KGML + model_name + ".xml"):
            copyfile(
                KGML + model_name + ".xml",
                KGML + model_name + "/" + model_name + ".xml",
            )

    def _create_dir_for_result(self):
        Path(RESULT).mkdir(parents=True, exist_ok=True)

    def _preprocessing(self, model):
        """
        start preprocessing from efmlrs
        """
        self._get_sbml(model)

        pre.start(
            KGML + model_name + "/" + model_name + ".sbml.xml",
            IGNORE_COMPARTIMENT,
            BOUNDS,
        )

    def _postprocessing_mplrs(self):
        """
        postprocessing mplrs
        """

        post.start(
            KGML + model_name + "/" + efms_mplrs,
            KGML + model_name + "/" + output_mplrs,
            KGML + model_name + "/" + info,
            False,
        )

    def _get_efms_from_file(self):
        efms = numpy.loadtxt(KGML + model_name + "/" + output_mplrs)
        efm_supps = [tuple(self.supp(efm)) for efm in efms]

        unique_efm_supps = set(efm_supps)
        self._efms = len(unique_efm_supps)

        return efms

    def _runKeggTranslatorToSbml(self):
        """
        KEGGtranslator can be used to batch convert multiple KGML formatted XML-files.
        To use this feature, simply give a folder as input argument, together with the desired output format and optional additional options.
        Example: java -jar KEGGtranslator.jar --input C:/ --format SBML
        """
        subprocess.check_output(
            [
                "java",
                "-jar",
                KEGGtranslator,
                "--input",
                KGML,
                "--format",
                "SBML_L3V1",
            ]
        )

    def _getSpecificKeggPathway(self, org="hsa"):
        """
        lists pathway, homo sapiens is default
        """
        api_url = "http://rest.kegg.jp/list/pathway/" + org
        response = requests.get(api_url)
        response.raise_for_status()
        return self._getOrganisms(response.text, org)

    def _getAllKeggOrganism(self):
        """
        get all available organisms from kegg
        """
        api_url = "http://rest.kegg.jp/list/organism/"
        response = requests.get(api_url)
        response.raise_for_status()
        # self._getOrganisms(response.text)

        org_split = response.text.split("\n")
        org = [i.split("\t")[1] for i in org_split if i]
        return org

    def _getOrganisms(self, response, org):
        """
        returns list of organisms
        """
        organisms = []

        r_split = response.split()
        for entry in r_split:
            if entry.startswith(org):
                # org = entry.split(":")[1]
                organisms.append(entry)
        return organisms

    def _getKGML(self, model):
        """
        get kgml file from KEGG
        """
        api_url = "http://rest.kegg.jp/get/" + model + "/kgml"
        file = KGML + model + ".xml"
        if not os.path.isfile(file):
            response = requests.get(api_url)
            response.raise_for_status()

            if not os.path.exists("KGML"):
                os.mkdir("KGML")
            with open(KGML + model + ".xml", "w") as f:
                f.write(response.text)

    def _getAllKGMLFromPathways(self, pathways):
        """
        From KEGG database
        Download all pathways for one organims (HSA default) as KGML file
        """
        for pathway in pathways:
            api_url = "http://rest.kegg.jp/get/" + pathway + "/kgml"
            file = KGML + pathway + ".xml"
            if not os.path.isfile(file):
                print("Downloading ", pathway)
                response = requests.get(api_url)
                response.raise_for_status()

                if not os.path.exists("KGML"):
                    os.mkdir("KGML")
                with open(KGML + pathway + ".xml", "w") as f:
                    f.write(response.text)
            else:
                print("Skipping ", pathway, " download. File already exists.")


parser = ArgumentParser()
parser.add_argument("-f", "--file", help="path to sbml file", metavar="SBML FILE")
parser.add_argument("-m", "--method", help="choose efmtool or efmlrs", metavar="METHOD")
parser.add_argument(
    "-p",
    "--pathway",
    help="download all pathway models for an organism from KEGG database, e.g. hsa for Homo sapiens (human)",
    metavar="ORGANISM",
)
parser.add_argument(
    "-k", "--keggtranslator", action="store_true", help="use KEGGTranslator"
)
parser.add_argument(
    "-sm",
    "--singlemodel",
    metavar="MODEL NAME",
    help="download single model from KEGG database",
)

args = parser.parse_args()


if __name__ == "__main__":
    e = efms()
    e._safety_checks()

    if args.singlemodel:
        e._getKGML(args.singlemodel)

    if args.pathway:
        pathway_res = e._getSpecificKeggPathway(args.pathway)
        e._getAllKGMLFromPathways(pathway_res)

    if args.keggtranslator:
        e._runKeggTranslatorToSbml()

    if args.file:
        # e.g. from file eco00030.sbml.xml, model_name will be eco00030,
        #     from file eco00030.xml, model_name will also be eco00030
        filename = os.path.basename(args.file)
        filename2 = os.path.splitext(filename)
        model_name = os.path.splitext(filename2[0])[0]
        MODEL = args.file

        if args.file.endswith(".sbml.xml"):
            efms_efmtool = model_name + ".sbml_cmp_efmtool.efms"
            efms_mplrs = model_name + ".sbml_cmp_mplrs.efms"
            info = model_name + ".sbml.info"
            output_mplrs = model_name + ".sbml_decomp_mplrs.efms"
            ine_cmp = KGML + model_name + "/" + model_name + ".sbml_cmp.ine"

        else:
            efms_efmtool = model_name + "_cmp_efmtool.efms"
            efms_mplrs = model_name + "_cmp_mplrs.efms"
            info = model_name + ".info"
            output_mplrs = model_name + "_decomp_mplrs.efms"
            ine_cmp = KGML + model_name + "/" + model_name + "_cmp.ine"

    if args.file and args.method is None:
        parser.error("--file requires --method.")

    if args.file and args.method == "efmtool":
        
        e._checkForReactions(MODEL)

        e._create_dir_for_model()
        h = hpy()

        t1_start = time()
        e._calc_efms_efmtool(MODEL)
        t1_stop = time()
        efmtool_totaltime = t1_stop - t1_start

        heap = h.heap()
        print("Elapsed time during calc with efmtool in seconds:", efmtool_totaltime)

        print("Total memory size: ", heap.size)

        result = {
            "model": model_name,
            "method": args.method,
            "metabolites": e.get_metabolites(),
            "reactions": e.get_reactions(),
            "reversibilities": e.get_reversibilities(),
            "irreversibilities": e.get_irreversibilities(),
            "efms": e.get_efms(),
            "runtime": efmtool_totaltime,
            "memory": heap.size,
        }

        e._create_dir_for_result()
        with open(RESULT + "efmtool_result_" + model_name + ".json", "w") as f:
            json.dump(result, f, indent=4)

    if args.file and args.method == "efmlrs":
        e._checkForReactions(MODEL)
        e._create_dir_for_model()
        h = hpy()

        pre1_start = time()
        e._preprocessing(MODEL)
        pre1_stop = time()
        pre_totaltime = pre1_stop - pre1_start

        print("Elapsed time during the preprocessing in seconds:", pre_totaltime)

        t1_start = time()
        e._calc_efms_mplrs(LRS, ine_cmp, KGML + model_name + "/" + efms_mplrs)
        t1_stop = time()
        mplrs_totaltime = t1_stop - t1_start

        print("Elapsed time during calc efms with mplrs in seconds:", mplrs_totaltime)

        post1_start = time()
        e._postprocessing_mplrs()
        post1_stop = time()
        post_totaltime = post1_stop - post1_start

        print("Elapsed time during the postprocessing in seconds:", post_totaltime)

        mpirun_totaltime = pre_totaltime + mplrs_totaltime + post_totaltime
        heap = h.heap()
        print("Total memory size: ", heap.size)

        e._get_efms_from_file()
        result = {
            "model": model_name,
            "method": args.method,
            "metabolites": e.get_metabolites(),
            "reactions": e.get_reactions(),
            "reversibilities": e.get_reversibilities(),
            "irreversibilities": e.get_irreversibilities(),
            "efms": e.get_efms(),
            "pre_runtime": pre_totaltime,
            "mplrs_totaltime": mplrs_totaltime,
            "post_totaltime": post_totaltime,
            "total_runtime": mpirun_totaltime,
            "memory": heap.size,
        }

        e._create_dir_for_result()
        with open(RESULT + "mplrs_result_" + model_name + ".json", "w") as f:
            json.dump(result, f, indent=4)
