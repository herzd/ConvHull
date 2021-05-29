
import chm_exact

reactions = [0, 1]

data_path = "../../DATA/toy/"

chm_exact.compute_CH(data_path + "toy_model_r.txt", data_path + "toy_model_S.txt", data_path + "toy_model_d.txt", reactions)