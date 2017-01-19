from aids import AID
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time


def get_aid_name(aid_string,f_name):
    aid_x = aid_string.split()
    aid_list = []
    name_list = []
    for bid in aid_x:
        aid_list.append(bid.replace("X",""))
    for cid in aid_list:
        time.sleep(6)
        aid = AID(cid)
        time.sleep(6)
        print(cid)
        name_list.append(''.join(aid.get_property("Name")))
        print(name_list)
    df = pd.DataFrame({"aid": aid_list,"name":name_list})
    print(df)
    df.to_csv(f_name + '.csv')

db15_441 = "X504749	X504834	X2701	X2707	X1458	X1490	X155	X157	X161	X165	X167	X175	X894	X884	X1460	X2546	X2551	X624044	X624297	X624030	X1053175	X902	X1478	X1030	X1688	X624032	X624296	X995	X504847"
#get_aid_name(db15_441,"db15_441")
db15_4 = "X115	X119	X123	X125	X7	X111	X117	X139	X15	X63	X79	X99	X113	X121	X27	X37	X41	X61	X65	X67	X73	X101	X105	X155	X157	X161	X165	X167	X175	X127	X409954	X409956	X17	X23	X81	X83"
get_aid_name(db15_4,"db15_4")