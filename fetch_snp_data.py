import requests,json,re
import threading
ENSEMBL_REST = "https://rest.ensembl.org"

class Fetch_Data(threading.Thread):
    def __init__(self,rsids,flank_length=800):
        super().__init__()
        self.result={}
        self.flank_length=flank_length
        self.rsids=rsids

    def run(self):
        r=requests.post(ENSEMBL_REST+"/variation/homo_sapiens", headers={ "Content-Type" : "application/json", "Accept" : "application/json"}, data=json.dumps({"ids":self.rsids}))
        r.raise_for_status()
        info=r.json()
        for rsid,v in info.items():
            possible_alleles=[i for i in v["mappings"][0]["allele_string"].strip('-/').split("/") if len(i)==1]
            if len(possible_alleles)<2:
                print(rsid,"has no valid alleles/less than 2 to choose from")
                continue
            p=int(v["mappings"][0]["start"])
            self.result[rsid]={
                "rsid":rsid,
                "chom":v["mappings"][0]["seq_region_name"],
                "allele_str":possible_alleles,
                "start":(t:=max(1,p-800)),
                "end":p+800,
                "rel_pos":self.flank_length if t>1 else p
                } 

def get_snp_data(rsid_info,lock,snp_list,info):
    seq_url = f"{ENSEMBL_REST}/sequence/region/human"
    seq_resp = requests.post(seq_url, headers={ "Content-Type" : "application/json", "Accept" : "application/json"},data=json.dumps({"regions" :rsid_info}))
    seq_resp.raise_for_status()
    resp = seq_resp.json()
    for i in resp:
        tre=info[i["query"]]
        for allele in tre["allele_str"]:
            # Validate that allele contains valid DNA characters only
            if not re.fullmatch("[ACGTNacgtn]+", allele):
                print(f"Skipping non-standard allele '{allele}' for {tre["rsid"]}")
                continue
            modified_seq = i["seq"][:tre["rel_pos"]] + allele.upper() + i["seq"][tre["rel_pos"] + 1:]
            if modified_seq[0]!='N':
                with lock:
                    snp_list.append({
                        "snpID": tre["rsid"],
                        "allele": allele.upper(),
                        "sequence": modified_seq,
                        "position": tre["rel_pos"]
                    })


def get_data(rsids,flank_length=800):
    all={}     
    threads=[Fetch_Data(rsids[i1:i1+200],flank_length) for i1 in range(0,len(rsids),200)]
    for i2 in threads:
        i2.start()   
    for i3 in threads:
        i3.join()
        all|=i3.result
    lock=threading.Lock()
    snp_list=[]
    info={}
    l1=[]
    for i in all.values():
        l1.append(t:=f"{i["chom"]}:{i["start"]}..{i["end"]}:1")
        info[t]=i
    if len(l1)>50:
        l2=[threading.Thread(target=get_snp_data,args=(l1[i:i+50],lock,snp_list,info)) for i in range(0,len(l1),50)]
        for i in l2:
            i.start()
        for i in l2:
            i.join()
    else:
        get_snp_data(l1,lock,snp_list,info)
    return snp_list, all


