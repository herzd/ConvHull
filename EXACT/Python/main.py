
import chm_exact

reactions = [0, 1]

data_path = "../../DATA/toy/"

chm_exact.compute_CH(data_path + "toy_reactions.txt", data_path + "toy_stoichs.txt", data_path + "toy_domains.txt", reactions)