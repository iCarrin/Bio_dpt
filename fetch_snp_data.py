import requests,json,re
from concurrent.futures import ThreadPoolExecutor,as_completed
ENSEMBL_REST = "https://rest.ensembl.org"

# https://rest.ensembl.org/variation/homo_sapiens/rs1799971

def request_data(url,data):
    r=requests.post(url, headers={ "Content-Type" : "application/json", "Accept" : "application/json"}, data=data)
    r.raise_for_status()
    return r.json()

def process_data(info,flank_length):
    for rsid,v in info.items():
        possible_alleles=[i for i in v["mappings"][0]["allele_string"].strip('-/').split("/") if len(i)==1]
        if len(possible_alleles)<2:
            print(rsid,"has no valid alleles/less than 2 to choose from")
            continue
        p=int(v["mappings"][0]["start"])
        yield {
            "rsid":rsid,
            "chom":v["mappings"][0]["seq_region_name"],
            "allele_str":possible_alleles,
            "start":(t:=max(1,p-800)),
            "end":p+800,
            "rel_pos":flank_length if t>1 else p
            } 
        
def get_snp_data(resp,info):
    for i in resp:
        tre=info[i["query"]]
        for allele in tre["allele_str"]:
            # Validate that allele contains valid DNA characters only
            if not re.fullmatch("[ACGTNacgtn]+", allele):
                print(f"Skipping non-standard allele '{allele}' for {tre["rsid"]}")
                continue
            modified_seq = i["seq"][:tre["rel_pos"]] + allele.upper() + i["seq"][tre["rel_pos"] + 1:]
            if modified_seq[0]!='N':
                yield {
                        "snpID": tre["rsid"],
                        "allele": allele.upper(),
                        "sequence": modified_seq,
                        "position": tre["rel_pos"]
                    }


def get_data(rsids,flank_length=800):
    all={}     
    info={}
    l1=[]
    with ThreadPoolExecutor() as e:
        f=[e.submit(request_data,f"{ENSEMBL_REST}/variation/homo_sapiens",json.dumps({"ids":rsids[i1:i1+200]})) for i1 in range(0,len(rsids),200)]
        for fe in as_completed(f):
            for i in process_data(fe.result(),flank_length):
                all[i["rsid"]]=i
                l1.append(t:=f"{i["chom"]}:{i["start"]}..{i["end"]}:1")
                info[t]=i

        snp_list=[]
        r=[e.submit(request_data,f"{ENSEMBL_REST}/sequence/region/human",json.dumps({"regions" :l1[i:i+50]})) for i in range(0,len(l1),50)]
        for re in as_completed(r):
            snp_list.extend(get_snp_data(re.result(),info))
                
    return snp_list, all


