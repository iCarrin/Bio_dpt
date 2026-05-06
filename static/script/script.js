const name = document.location.host
console.log(name)
let url;
document.getElementById("submit").addEventListener('click', async () => {
    if (!document.getElementById("file-upload").files[0] && !document.getElementById("input-field").value){
        alert("no values inputed or file uploaded")
        return
    }
    obj = {}
    l1 = ["mv_conc", "dv_conc", "dntp_conc", "dna_conc", "dmso_conc", "dmso_fact", "formamide_conc", "annealing_temp_c", "temp_c", "tm_method", "salt_corrections_method", "desired_tm", "diff", "homodimer_goal", "hairpin_goal", "target_gc", "heterodimer_max", "primer_dimer_distance", "min_probe_len", "max_probe_len", "min_primer_len", "max_primer_len", "min_primer_dist", "max_primer_dist", "flank_length", "pdfoutput_precision"]
    l1.forEach(id => {
        obj[id] = document.getElementById(id.replaceAll("_", "-")).value

    });
    console.log(obj)
    const formData = new FormData()
    formData.append("file", document.getElementById("file-upload").files[0])
    formData.append("input", document.getElementById("input-field").value)
    formData.append("settings", JSON.stringify(obj))

    responce = await fetch(`http://${name}/api/gids`, {
        method: 'POST', body: formData
    })

    if (responce.ok) {
        t1=await responce.blob()
        url=URL.createObjectURL(t1)
        document.getElementById("pdfFrame").src = url;    
    }

})
