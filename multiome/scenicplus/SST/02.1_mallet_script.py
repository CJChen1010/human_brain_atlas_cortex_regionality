from pycisTopic.lda_models import *
import os
import sys
os.environ["MALLET_MEMORY"] = "250G"

def main(n_topic, obj_path):
    cistopic_obj = pickle.load(open(obj_path, "rb"))
    out_dir = "/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/SST/pycistopic_outputs"
    path_to_mallet_binary = "/tscc/nfs/home/biy022/softwares/scenicplus/Mallet-202108/bin/mallet"
    tmp_path = "/tscc/projects/ps-epigen/users/biy022/biccn/data/scplus_tmp/{}".format(n_topic)
    os.makedirs(tmp_path, exist_ok=True)
    model = run_cgs_models_mallet(
        mallet_path=path_to_mallet_binary,
        cistopic_obj=cistopic_obj,
        n_topics=[int(n_topic)],
        n_cpu=20,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=tmp_path,
        save_path=out_dir
    )

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
