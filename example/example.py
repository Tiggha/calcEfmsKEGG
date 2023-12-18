from EFM_AIO import *

KGML = "models/"

if __name__ == "__main__":
    e = efms()

    if args.singlemodel:
        e._getKGML(args.singlemodel)

    if args.pathway:
        pathway_res = e._getSpecificKeggPathway(args.pathway)
        e._getAllKGMLFromPathways(pathway_res)

    if args.keggtranslator:
        e._runKeggTranslatorToSbml()

    if args.file:
        # e.g. from file eco00030.sbml.xml, MODEL_NAME will be eco00030,
        #     from file eco00030.xml, MODEL_NAME will also be eco00030
        filename = os.path.basename(args.file)
        filename2 = os.path.splitext(filename)
        MODEL_NAME = os.path.splitext(filename2[0])[0]
        MODEL = args.file

        if args.file.endswith(".sbml.xml"):
            e.set_fileNames(True, MODEL_NAME)
            EFMS_MPLRS = MODEL_NAME + ".sbml_cmp_mplrs.efms"
            INE_CMP = KGML + MODEL_NAME + "/" + MODEL_NAME + ".sbml_cmp.ine"

        else:
            e.set_fileNames(False, MODEL_NAME)
            EFMS_MPLRS = MODEL_NAME + "_cmp_mplrs.efms"
            INE_CMP = KGML + MODEL_NAME + "/" + MODEL_NAME + "_cmp.ine"

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
            "model": MODEL_NAME,
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
        with open(RESULT + "efmtool_result_" + MODEL_NAME + ".json", "w") as f:
            json.dump(result, f, indent=4)

    if args.file and args.method == "efmlrs":
        e._checkForReactions(MODEL)
        e._create_dir_for_model(MODEL_NAME)
        h = hpy()

        pre1_start = time()
        e._preprocessing(MODEL)
        pre1_stop = time()
        pre_totaltime = pre1_stop - pre1_start

        print("Elapsed time during the preprocessing in seconds:", pre_totaltime)

        t1_start = time()
        e._calc_efms_mplrs(LRS, INE_CMP, KGML + MODEL_NAME + "/" + EFMS_MPLRS)
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
            "model": MODEL_NAME,
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
        with open(RESULT + "mplrs_result_" + MODEL_NAME + ".json", "w") as f:
            json.dump(result, f, indent=4)
