import sys
import os
from scipy.spatial import distance
from pymongo import MongoClient
import json
# 使用过滤过的结果做相似性柱状图
import pyecharts
from pyecharts import Page
from pyecharts import Bar
#import pyecharts.echarts.events as events
#from pyecharts_javascripthon.dom import alert

incoming_smiles = sys.argv[1] # read incoming molecule from console
query_db = sys.argv[2]
similarity_method = sys.argv[3] # read incoming similarity method
similarity_min = float(sys.argv[4])
similarity_max = float(sys.argv[5])
jobID = sys.argv[6]
username = sys.argv[7]


# 结果保存路径/home/zyli/webserver/runjobs/find_similar_mols/username/jobid
save_path = "/home/zyli/webserver/runjobs/find_similar_mols/" + username + "/" + jobID + "/python_results/"
if not os.path.exists(save_path):
    os.makedirs(save_path)

try:
    #conn = MongoClient(host='localhost', port=27018)  # init a connection
    conn = MongoClient(host='localhost', port=27018, username='brs3d', password='j0nathan!', authSource='admin')
    my_db = conn.local  # connnect mongodb, database: local

    # query molecule zscore
    ref_raw_result = my_db.allZScore.find_one({"smiles" : incoming_smiles})

    # connnect collection brsData
    if query_db == "chemdiv":
        my_collection = my_db.chemdivZScore
    elif query_db == "enamine":
        my_collection = my_db.enamineZScore
    elif query_db == "zinc":
        my_collection = my_db.zincZScore
    elif query_db == "gdd":
        my_collection = my_db.gddZScore
    elif query_db =="gll":
        my_collection = my_db.gllZScore
    else:
        my_collection = my_db.specsZScore

    if ref_raw_result == None:
        # send submitting online calculation signal
        print("no matching result found in the database")

    else:
        all_result_dict = {}
        the_rest_of_db_raw_results = my_collection.find({'smiles' : {'$ne' : incoming_smiles}})#.limit(10000)
        #print(the_rest_of_db_raw_results)
        if similarity_method == 'cosine':
            for each_raw_result_dic in the_rest_of_db_raw_results:
                compare_set ={"ref_raw_result":ref_raw_result, "each_raw_result_dic": each_raw_result_dic}

                list1 = list(compare_set['ref_raw_result'].values())
                list2 = list(compare_set['each_raw_result_dic'].values())
                target_smiles = list(compare_set['each_raw_result_dic'].values())[-1]
                dist = distance.cosine(list1[2:301], list2[2:301])
                sim = 1 - dist
                all_result_dict[target_smiles] = sim

        if similarity_method == 'correlation':
            for each_raw_result_dic in the_rest_of_db_raw_results:
                compare_set ={"ref_raw_result":ref_raw_result, "each_raw_result_dic": each_raw_result_dic}

                list1 = list(compare_set['ref_raw_result'].values())
                list2 = list(compare_set['each_raw_result_dic'].values())
                target_smiles = list(compare_set['each_raw_result_dic'].values())[-1]
                dist = distance.correlation(list1[2:301], list2[2:301])
                sim = 1 - dist
                all_result_dict[target_smiles] = sim

        # 按照相似性值降序排列
        sorted_all_result_list_by_value = sorted(all_result_dict.items(), key=lambda d: d[1])
        #print("%s" %sorted_all_result_list_by_value, file=open(save_path + query_db+ "_all_sim_result.txt", "w"))

        # 按照相似性阈值筛选数组中的结果
        all_result_count = len(sorted_all_result_list_by_value)
        # 使用列表推导式过滤元组中的相似性值
        within_threshold_result = [item for item in sorted_all_result_list_by_value if
                                   (item[1] > (similarity_min/100) and item[1] <= (similarity_max/100))]
        # 统计阈值内的分子个数
        within_threshold_count = len(within_threshold_result)

        # 将过滤的结果转化为json格式，提供给node.js中的controller使用
        within_threshold_result_json = json.dumps(within_threshold_result)
        print("%s" %within_threshold_result_json, file=open(save_path +query_db+ "_within_threshold_sim.json", "w"))

        # todo:将整个数据库相似性总结果也转化为json保存到硬盘
        sorted_all_result_list_by_value_json = json.dumps(sorted_all_result_list_by_value)
        print("%s" %sorted_all_result_list_by_value_json, file=open(save_path + query_db+ "_sorted_all_result_list_by_value.json", "w"))

        # 使用pyecharts作图
        all_compare_smiles = []
        all_compare_sim = []
        for i in within_threshold_result:
            compare_smiles = i[0]
            compare_sim = i[1]
            all_compare_smiles.append(compare_smiles)
            all_compare_sim.append(compare_sim)
        # 横轴
        attr = [i for i in all_compare_smiles]
        # 纵轴
        v1 = [j for j in all_compare_sim]
        
        bar = Bar("BRS-3D - find similar mols plot", 'Total '+ str(within_threshold_count) +' similar molecules found')
        bar.add("", attr, v1, yaxis_name='similarity', yaxis_name_gap=38, bar_category_gap='20%', label_color=["#4e79a7"], is_label_show=False, is_xaxis_show=False, is_datazoom_show=True, is_stack=False, is_more_utils=True)
        # 定义渲染出来的webpage路径
        render_path = save_path + query_db+ '_brs3d_sim_plot.html'
        bar.render(render_path)


except Exception as e:
    print(e)
