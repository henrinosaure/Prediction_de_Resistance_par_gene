# Author : E. Fortin
# Date : month 20--
# 
# 
#    Description
#  ═══════════════
#  - One_line_description -
# 
# Further_description
# 
# 
# 
#    Arguments
#  ═════════════
# • input_dir (String)
#     Absolute or relative path of input directory containing input files.
# 
# • output_dir (String)
#     Absolute or relative path of output directory.
# 
# • -o | --opt_arg (String) [default: default_value]
#     The_optional_argument
# 
# 
#    Outputs
#  ═══════════
# bla_bla
# 
# 
#    Returns
#  ═══════════
# bla_bla
#
# %%════════════════════════════════════════╗
# ║                  Main                   ║
# ╚═════════════════════════════════════════╝

plier <- function(p_len_data, p_nplis, p_val_ratio) {
    idx = sample(1:p_len_data)
    
    groupes <- split(idx, cut(seq_along(idx), p_nplis, labels = FALSE))
    
    plis = vector("list", p_nplis)
    
    for(i in 1:p_nplis) {
        
        groupe <- groupes[[i]]
        nb_obs <- length(groupe)
        nb_val <- round(nb_obs * p_val_ratio)
        
        plis[[i]]$val <- groupe[1 : nb_val]
        plis[[i]]$train <- groupe[(nb_val+1) : nb_obs]
    }
    
    plis
}