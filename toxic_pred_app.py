
# libraries --------------------------------------------------------
import pickle
from numpy.core.fromnumeric import shape
import streamlit as st
import SessionState
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
from datetime import date
import dash
import dash_table
import dash_core_components as dcc

# libraries for pfeature-feature calculations 
from Pfeature.pfeature import *
from Pfeature.pfeature import aac_wp
import sys
import importlib
importlib.reload(sys)
if sys.version[0] == '2':
    importlib.reload(sys)
    sys.setdefaultencoding("utf-8")
from six import StringIO

# libraries for xgboost 
from sklearn.metrics import classification_report, accuracy_score
import xgboost
from sklearn.metrics import accuracy_score
from xgboost import XGBClassifier
from sklearn.datasets import load_iris
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, KFold
from sklearn import metrics
from  sklearn.datasets import make_hastie_10_2
from xgboost.sklearn import XGBClassifier


# Set to use full page ----------------------------------------------
st. set_page_config(layout="wide")

#title --------------------------------------------------------------
st.title('Peptides Toxicity Prediction  -   多肽毒性预测')

st.write('\n')
st.write('\n')

# Intro ------------------------------------------------------------
st.markdown("<br>", unsafe_allow_html=True)
"""使用步骤:\n
1. 您可以选择直接粘贴fasta格式的多肽序列至下方文本框，或上传一个包含多肽序列的fasta文档。请注意，一次只能使用一种方式！ \n
2. 如果输入多肽可以被计算，“您的多肽特征已经计算完毕 ，准备开始预测毒性” 的指示将会出现。
3. 多肽毒性预测的结果将会被显示在最下方，毒性越大，数字越大。如需下载结果，请点击最后的下载按钮。
---
"""

# Global Variable ---------------------------------------------------
store_tol = []
store_text = []
store_sequence_in = []
process_text = '多肽特征正在计算中...'
input_name = ''


# Method 1: Input sequence -------------------------------------------

# An example
default_input_sequence = ">test seq 1 （有毒）|\nEDGYLLNRDTGCKVSCGTCRYCND"
user_input = st.text_area("◉ 方式1：请在下面的方框输入您的fasta格式多肽序列, 完成后请按 ctrl+C 或用鼠标点击屏幕空白处保存 （下方是一个案例序列）", default_input_sequence)


# file location
fw = open("save_input_seq.fa", 'w')
split_input = user_input.split("\n") 
for line in split_input:  
    #write strings for views  
    fw.write(line.rstrip("\n"))    
    fw.write("\n") 

fw.close()

#!!!!!!!!!!!!!!!!!!!!!目前还未解决 > 会变成自动格式的问题
for line in split_input:  
    #write strings for views  
    line_now = line.rstrip("\n")
    store_tol.append(line_now)
    #st.write(line_now)    

input_name = "save_input_seq.fa"

st.write('\n')

# Method 1: Input a file directly --------------------------------------------------
uploaded_file = st.file_uploader("◉ 方式2：请选择打开fasta格式（.fa）的文件夹/文件")
if uploaded_file is not None:
    # To read file as bytes:
    bytes_data = uploaded_file.getvalue()
    # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    # To read file as string:
    string_data = stringio.read()

    if uploaded_file is not None:
        store_tol = []
        store_text = []
        store_sequence_in = []
        # input_name = input_filename
        # file1 = open(input_name, 'r')
        # Lines = file1.readlines()

        fw = open("save_input_file.fa", 'w')

        #st.write(type(split_input)) 
        split_input = string_data.split("\n")  
        for line in split_input:  
            #write strings for views  

           line_now = line.rstrip("\n\t")
           if line_now[-1].isalpha() is False: 
               line_now = line_now[:-1]   
           fw.write(line_now) 
           fw.write("\n") 

        fw.close()
        input_name = "save_input_file.fa"
     
        for line in split_input:  
        #write strings for views  
            line_now = line.rstrip("\n\t")
            if line_now[-1].isalpha() is False: 
                line_now = line_now[:-1]
            store_tol.append(line_now)
        
    else:
        input_name = "save_input_seq.fa"

st.write('\n')
st.write('\n')
st.write('您的输入为：\n')

for line in split_input:  
    #write strings for views  
    line_now = line.rstrip("\n")
    #store_tol.append(line_now)
    st.write(line_now)   

#文件必须要在跑这个streamlit app的同文件夹更深处
# def file_selector(folder_path='.', target="文件夹/文件"):
#     filenames = [f for f in os.listdir(folder_path) if
#                  not f[0] == "."]  # get file names from dir excluding hidden files
#     selected_filename = st.selectbox(f'方式2：请选择打开fasta格式（.fa）的{target}', filenames)
#     abs_path = os.path.join(folder_path, selected_filename)
#     if os.path.isdir(abs_path):
#         return file_selector(abs_path, target)
#     return os.path.join(folder_path, selected_filename)

# input_filename = file_selector()
# st.write('已选择 `%s`' % input_filename)

# button widget ---------------------------------------------------

st.write('\n')
st.write('\n')
session_state = SessionState.get(checkboxed=False)
if st.button('点击开始预测毒性') or session_state.checkboxed:
    session_state.checkboxed = True
    st.write('\n')
    st.write('\n')
    st.write(process_text)

    store_text = store_tol[0::2]
    store_sequence_in = store_tol[1::2]

    # Process data -------------------------------------------------------

    # Amino-acid Composition
    output_name_1 = "save_input_seq_feature_AAC.csv"
    aac_wp(input_name, output_name_1)

    # Amino-acid Composition of C-Terminal
    output_name_21 = "save_input_seq_feature_C5AAC.csv"
    aac_ct(input_name, output_name_21, 5)

    # Amino-acid Composition of N-Terminal
    output_name_22 = "save_input_seq_feature_N5AAC.csv"
    aac_nt(input_name, output_name_22, 5)

    # Dipeptide Composition
    output_name_3 = "save_input_seq_feature_DPC.csv"
    dpc_wp(input_name, output_name_3, 1)

    # Tripeptide Composition
    # output_name_4 = "save_input_seq_feature_TPC.csv"
    # tpc_wp(input_name, output_name_4)

    # Atom Composition
    output_name_5 = "save_input_seq_feature_ATC.csv"
    atc_wp(input_name, output_name_5)

    #read file
    f_1 = pd.read_csv(output_name_1) 
    f_21 = pd.read_csv(output_name_21) 
    f_22 = pd.read_csv(output_name_22) 
    f_3 = pd.read_csv(output_name_3) 
    #f_4 = pd.read_csv(output_name_4) 
    f_5 = pd.read_csv(output_name_5) 

    #concatenated dataframes
    conc_df_1 = pd.concat([f_1, f_21], axis=1)
    conc_df_2 = pd.concat([f_22, f_3], axis=1)
    #conc_df_3 = pd.concat([f_4, f_5], axis=1)

    tot_conc = pd.concat([conc_df_1, conc_df_2], axis=1)
    tot_conc = pd.concat([tot_conc, f_5], axis=1)

    #save file
    tot_conc.to_csv('total_input_features.csv', header=True, index=False)

    # Start Predictions----------------------------------------------
    process_text = '您的多肽特征已经计算完毕 ，准备开始预测毒性... \n'
    st.write(process_text)

    # load saved model
    filename = 'finalized_model.sav'
    clf = pickle.load(open(filename, 'rb'))

    # load the features of the seqs that need to be predicted 
    X_test = pd.read_csv("total_input_features.csv")
    # get prediction result
    size_result = X_test.size / 465
    y_pred=clf.predict_proba(X_test).T[1]
    y_pred = pd.to_numeric(y_pred)

    st.write('您的多肽毒性已经预测完毕 ，结果如下： \n')

    y_print_id = y_pred
    y_id = []
    i = 0

    # If need the prediction result to show the txt
    # instead of the probabilty instead
    # while i < size_result:
    #     now_num = y_pred[i].item()
    #     # st.write(type(now_num))
    #     y_print_id[i] = i + 1
    #     if now_num == 0:
    #         y_print_result.append("无毒性") 
    #         # y_print_result[i][1] = "无毒性"
    #     else:
    #         y_print_result.append("有毒性")
    #         # y_print_result[i][1] = "有毒性"
    #     i = i + 1

    # Unify the type of all data
    while i < size_result:
        y_id.append(i + 1)
        now_num = y_pred[i].item()
        y_pred[i] = float("{0:.5f}".format(now_num)) 
        i = i + 1

    store_text = np.array(store_text)
    store_sequence_in = np.array(store_sequence_in)
    # Concatenate data
    test = np.stack((y_id, store_text, store_sequence_in, y_pred))
    turn_test = test.T
    # Store numpy data as dataframe
    # (For future file-saving purposes)
    store_as_df_0 = pd.DataFrame(turn_test)
    store_as_df_0.columns =['序号', '多肽序号', '多肽序列', '预测结果']

    # Plot the result -----------------------------------------------
    import plotly.graph_objects as go

    fig = go.Figure(data=[go.Table(
        columnorder = [1,2,3,4],
        columnwidth = [60,150,400,140],
        header=dict(values=['<b>序号</b>', '<b>多肽序号</b>', '<b>多肽序列</b>', '<b>预测结果</b>']),
                    cells=dict(values=test))
                        ])

    st.write(fig)

    # Click button to save results as a csv file ----------------------
    # Static state was used to save the results when
    # clicking a button in button
    if st.button('点击下载多肽毒性预测结果（.csv）'):
        today = date.today()
        now = datetime.now()
        d4 = today.strftime("%m-%d-%Y")
        current_time = now.strftime("%H:%M:%S")
        dl_name = 'toxicity_pred_' + d4 + '_' + current_time + '.csv'
        store_as_df_0.to_csv(dl_name,index=False, encoding='utf_8_sig')
        string_print = '您的多肽毒性预测结果已经下载完毕 ，名字为： ' + dl_name
        st.write(string_print)


