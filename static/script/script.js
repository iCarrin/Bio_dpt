const name = document.location.host
console.log(name)
let url;

function create_spinner() {
    const spinner = document.createElement("div");

    // Style the spinner
    spinner.style.width = "50px";
    spinner.style.height = "50px";
    spinner.style.border = "6px solid #ddd";
    spinner.style.borderTop = "6px solid #3498db";
    spinner.style.borderRadius = "50%";
    spinner.style.animation = "spin 1s linear infinite";
    spinner.style.margin = "20px auto";
    spinner.id="spinner"

    // Add spinner to page
    // document.getElementById("output").appendChild(spinner);
    document.getElementById("output").prepend(spinner)

    // Create animation using JS
    const style = document.createElement("style");
    style.textContent = `
@keyframes spin {
  from {
    transform: rotate(0deg);
  }
  to {
    transform: rotate(360deg);
  }
}
`;

    document.head.appendChild(style);
}

function remove_spinner() {
     document.getElementById("spinner").remove()
}

document.getElementById("submit").addEventListener('click', async () => {
    create_spinner()
    if (!document.getElementById("file-upload").files[0] && !document.getElementById("input-field").value) {
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
    remove_spinner()
    if (responce.ok) {
        t1 = await responce.blob()
        url = URL.createObjectURL(t1)
        document.getElementById("pdfFrame").src = url;
    }

})
