from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    TableStyle,
    Table
)
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from Primer_Classes import Probe,Primer

def json_to_pdf(json_data, output_file,pagesize=(900,375)):
    doc = SimpleDocTemplate(output_file, pagesize=pagesize)
    styles = getSampleStyleSheet()
    story = []
    # Title
    if "title" in json_data:
        story.append(Paragraph(json_data["title"], styles["Title"]))
        story.append(Spacer(1, 20))
    # Elements
    for element in json_data.get("elements", []):
        if element["type"] == "paragraph":
            story.append(Paragraph(element["text"], styles["Normal"]))

        elif element["type"] == "spacer":
            story.append(Spacer(1, element.get("height", 12)))

        elif element["type"] == "table":
            table = Table(element["data"])
            table.setStyle(TableStyle([
                ("GRID", (0,0), (-1,-1), 1, colors.black),
                ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                ("ALIGN", (0,0), (-1,-1), "CENTER"),
            ]))
            story.append(table)
            story.append(Spacer(1, 20))

    doc.build(story)

def create_output_json(info,output_file,page_size,percision=2):
    if len(info)==0:return
    d1={"elements":[{"type":"table","data":[]}]}
    if type(info[0])==Probe or type(info[0])==Primer:
        d1["elements"][0]['data'].append([i for i in vars(info[0])])
    else:
        d1["elements"][0]['data'].append([i for i in info[0]])
        
    for i in info:
        if type(i)==Probe  or type(i)==Primer:
            d1["elements"][0]['data'].append(i.to_list(percision))
        else:
            d1["elements"][0]['data'].append(list(i.values()))
    json_to_pdf(d1,output_file,page_size)


        



