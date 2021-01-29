import bbh
# check readme file for more explanation
# create files with best hits against each other and against itself
CE_dbDP_besthists = bbh.retrieve_best_hist("queryCE_VS_databaseDP_ortho", "CEvsdbDP_besthits.txt")
DP_dbCE_besthists = bbh.retrieve_best_hist("queryDP_VS_databaseCE_ortho", "DPvsdbCE_besthits.txt")
CE_dbCE_besthists = bbh.retrieve_best_hist("query_paralogs_ce", "CE_dbCE_besthists.txt")
DP_dbDP_besthists = bbh.retrieve_best_hist("query_paralogs_dp", "DP_dbDP_besthists.txt")

# create file with orthologs and paralogs ; output file : Ce_Dp_orthologs.txt
bbh = bbh.bbh_2species("CEvsdbDP_besthits.txt", "DPvsdbCE_besthits.txt","CE_dbCE_besthists.txt", "DP_dbDP_besthists.txt")
