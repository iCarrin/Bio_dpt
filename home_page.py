from flask import Flask,Response, render_template,request,jsonify
import json
from fetch_snp_data import get_data
from Multiplexer_class import Multiplexer
l1=[
        ("mv_conc",float), ("dv_conc",int), 
        ("dntp_conc",float), ("dna_conc",int), 
        ("dmso_conc",float), ("dmso_fact",float), 
        ("formamide_conc",float), ("annealing_temp_c",float), 
        ("temp_c",float), ("tm_method",str),
        ("salt_corrections_method",str),
        ("desired_tm", float),("diff", float),
        ("homodimer_goal", float),("hairpin_goal", float),
        ("target_gc",float),
        ("heterodimer_max",float),("primer_dimer_distance",int),
        ("min_probe_len",int),("max_probe_len",int),
        ("min_primer_len",int),("max_primer_len",int),
        ("min_primer_dist",int),("max_primer_dist",int),
        ("flank_length",int),("pdfoutput_precision", int)
        ]

app = Flask(__name__)



@app.route("/")
def index():
    return render_template('index.html')

@app.route("/api/gids",methods=['POST'])
def get_ids():
    file=request.files.get('file')
    inp=request.form.get('input')
    settings=request.form.get('settings')
    d1=json.loads(settings)
    for k,t in l1:
        d1[k]=t(d1[k])
    s1=set()
    if file:
        for i in request.files.get('file').stream:
            s1.add(i.decode('utf-8').strip())
    if inp:
        for i in inp.split():
            s1.add(i.strip())
    
    snp_df, allele_df=get_data(rsids=list(s1),web=True)
    tre=Multiplexer(snp_df, allele_df,**d1).main(web=True)
    return Response(tre,mimetype="application/pdf",headers={"Content-Disposition": "inline"})
    


if __name__ == "__main__":
    app.run(port=8000)