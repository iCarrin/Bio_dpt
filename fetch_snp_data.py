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
        for k,v in info.items():
            l1=[i for i in v["mappings"][0]["allele_string"].strip('-/').split("/") if len(i)==1]
            if len(l1)<2:
                print(k,"has no valid alleles/less than 2 to choose from")
                continue
            p=int(v["mappings"][0]["start"])
            self.result[k]={
                "rsid":k,
                "chom":v["mappings"][0]["seq_region_name"],
                "allele_str":l1,
                "start":(t:=max(1,p-800)),
                "end":p+800,
                "rel_pos":self.flank_length if t>1 else p
                } 

def ask_user(rsid,raw_alleles):
        print(f"alleles for {rsid}:")
        for i, allele in enumerate(raw_alleles):
            print(f"{i+1}) {allele}")
        alleles_wanted = input("type the corresponding numbers for the alleles you want separated by spaces, or just \"All\" for all of them: ")
        if alleles_wanted.strip().upper() == "ALL" or alleles_wanted.strip().upper() == "":
            print("using all alleles")
            wanted_alleles= raw_alleles
            return wanted_alleles
        elif re.fullmatch(r'[0-9\s]+', alleles_wanted):
            indices = [int(x) for x in alleles_wanted.strip().split(" ")]
            if max(indices) > len(raw_alleles) or min(indices) < 1:
                print("you asked for an allele that wasn't in the list")
                return ask_user()
            else:
                wanted_alleles=[raw_alleles[i-1] for i in indices]
                return wanted_alleles
        else:
            print("you typed more than just numbers and spaces. Try again")
            return ask_user()


def get_snp_data(rsid_info,lock,snp_list,wanted_alleles):
    seq_url = f"{ENSEMBL_REST}/sequence/region/human/{rsid_info["chom"]}:{rsid_info["start"]}..{rsid_info["end"]}:1?"
    seq_resp = requests.get(seq_url, headers={"Content-Type": "text/plain"})
    seq_resp.raise_for_status()
    template_seq = seq_resp.text.strip()
    for allele in wanted_alleles:
            # Validate that allele contains valid DNA characters only
            if not re.fullmatch("[ACGTNacgtn]+", allele):
                print(f"Skipping non-standard allele '{allele}' for {rsid_info["rsid"]}")
                continue
            modified_seq = template_seq[:rsid_info["rel_pos"]] + allele.upper() + template_seq[rsid_info["rel_pos"] + 1:]
            if modified_seq[0]!='N':
                with lock:
                    snp_list.append({
                        "snpID": rsid_info["rsid"],
                        "allele": allele.upper(),
                        "sequence": modified_seq,
                        "position": rsid_info["rel_pos"]
                    })


def get_rsids():
    file=input("load from file (y/n) ")
    if file.lower() in ['y','yes']:
        filename= input("input file path to load ")
        with open(filename) as txt:
            return {i.strip() for i in txt}
    else:
        l1=set()
        while True:
            rsid=input("enter rsid type 'exit' to exit: ")
            if rsid.lower()!="exit":
                l1.add(rsid)
            else:
                break
        return l1


def get_data(flank_length=800):
    rsids=list(get_rsids())
    all={}     
    threads=[Fetch_Data(rsids[i:i+200],flank_length) for i in range(0,len(rsids),200)]
    for i in threads:
        i.start()
    for i in threads:
        i.join()
        all|=i.result
    lock=threading.Lock()
    snp_list=[]
    threads1=[(t:=threading.Thread(target=get_snp_data,args=(i,lock,snp_list,ask_user(i["rsid"],i["allele_str"]))),t.start())[0] for i in all.values()]    
    for i in threads1:
        i.join()
    return snp_list


   


