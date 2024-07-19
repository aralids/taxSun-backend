from flask import Flask, render_template, request, jsonify, json
import copy
import taxopy
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

@app.route('/')
def index():
  return render_template('index.html')

database_dir = "."
taxdb = taxopy.TaxDb()
rankPatternFull = ["root", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "superorder", "order", "suborder", "superfamily", "family", "subfamily", "supergenus", "genus", "subgenus", "superspecies", "species"]

@app.route('/load_tsv_data', methods=["POST"])
def load_tsv_data():
    if request.method == 'POST':
        print("request.files: ", request.files)
        print("request.data: ", request.data)
        file = request.files['file'].read()
        file_lines = (file.decode("utf-8")[:-1]).split("\n")
        header_line = file_lines[0]
        lines = file_lines[1:]
        
        raw_tax_set, raw_lns, e_value_enabled, fasta_enabled = calc_raw_tax_set(header_line, lines)
        tax_set, lns = calc_tax_set(raw_tax_set, raw_lns, e_value_enabled, fasta_enabled)
        lns = sort_n_uniquify(lns)
        tax_set = correct_tot_counts(lns, tax_set)
        print("len(lns): ", len(lns))

        return jsonify({"lns": lns, "taxSet": tax_set, "eValueEnabled": e_value_enabled, "fastaEnabled": fasta_enabled, "rankPatternFull": rankPatternFull})


def calc_raw_tax_set(header_line, lines):
    raw_tax_set = {"root root": {"taxID": "1", 
                                "rawCount": 0, 
                                "totCount": 0, 
                                #"lineageNames": [["root", "root"]], 
                                "name": "root",
                                "rank": "root",  
                                "lnIndex": 0,
                                "names": [["root root", -1]], 
                                "geneNames": [], 
                                "eValues": [], 
                                "fastaHeaders": [], 
                                "children": []}}
    id_set = {"1": "root root"}
    id_sum = 0
    raw_lns = [[["root", "root"]]]

    e_value_enabled = False
    if "evalue" in header_line.lower():
        e_value_enabled = True

    fasta_enabled = False
    if "fasta" in header_line.lower():
        fasta_enabled = True
    

    for line in lines:
        line2list = line.split("\t")
        line2list = [item.replace("\r","") for item in line2list]
        gene_name = line2list[0]
        taxID = line2list[1]

        if len(line2list) >= 3:
            e_value = line2list[2]

        if len(line2list) >= 4:
            fasta_header = line2list[3]

        if taxID == "NA" or taxID == "":
            taxID = "1"

        if not (taxID in id_set):
            id_sum += int(taxID)
            taxon = taxopy.Taxon(int(taxID), taxdb)
            name = taxon.name
            rank = taxon.rank
            lineageNamesList = taxon.rank_name_dictionary
            dictlist = [["root", "root"]] + [[k,v] for k,v in lineageNamesList.items()][::-1]
            if not (dictlist[-1][0] == rank and dictlist[-1][1] == name):
                dictlist += [[rank, name]]
            raw_tax_set[name + " " + rank] = {"taxID": taxID, 
                                              #"lineageNames": dictlist, 
                                              "rawCount": 1, 
                                              "name": name,
                                              "rank": rank, 
                                              "totCount": 1, 
                                              "names": [[name + " " + rank, 0]], 
                                              "geneNames": [gene_name], 
                                              "children": []}
            
            if len(line2list) >= 3:
                if e_value != "":
                    raw_tax_set[name + " " + rank]["eValues"] = [float(e_value)]
                else:
                    raw_tax_set[name + " " + rank]["eValues"] = [None]

            if len(line2list) >= 4:    
                if fasta_header != "":
                    raw_tax_set[name + " " + rank]["fastaHeaders"] = [fasta_header]
                else:
                    raw_tax_set[name + " " + rank]["fastaHeaders"] = [None]

            id_set[taxID] = name + " " + rank
            raw_lns.append(dictlist)
        else:
            raw_tax_set[id_set[taxID]]["totCount"] += 1
            raw_tax_set[id_set[taxID]]["rawCount"] += 1
            raw_tax_set[id_set[taxID]]["geneNames"].append(gene_name)
            raw_tax_set[id_set[taxID]]["names"][0][1] += 1

            if len(line2list) >= 3:
                if e_value != "":
                    raw_tax_set[id_set[taxID]]["eValues"].append(float(e_value))
                else:
                    raw_tax_set[id_set[taxID]]["eValues"].append(None)

            if len(line2list) >= 4:
                if fasta_header != "":
                    raw_tax_set[id_set[taxID]]["fastaHeaders"].append(fasta_header)
                else:
                    raw_tax_set[id_set[taxID]]["fastaHeaders"].append(None)
    return raw_tax_set, raw_lns, e_value_enabled, fasta_enabled

def calc_tax_set(raw_tax_set, raw_lns, e_value_enabled, fasta_enabled):
    existent = {}
    deleted ={}
    created = {}
    lns = copy.deepcopy(raw_lns)

    for i in reversed(range(0, len(raw_lns))):
        ln = raw_lns[i]

        inherited_taxon = ""
        inherited_unaCount = 0
        inherited_geneNames = []

        if e_value_enabled:
            inherited_eValues = []
        if fasta_enabled:
            inherited_fastaHeaders = []

        for j in reversed(range(0, len(ln))):
            name = ln[j][1]
            rank = ln[j][0]
            taxon = name + " " + rank

            if rank in rankPatternFull:
                if taxon in raw_tax_set:
                    if not (taxon in existent):
                        existent[taxon] = raw_tax_set[taxon]
                        existent[taxon]["unaCount"] = raw_tax_set[taxon]["rawCount"]
                    if inherited_unaCount > 0:
                        existent[taxon]["unaCount"] += inherited_unaCount
                        existent[taxon]["totCount"] += inherited_unaCount
                        existent[taxon]["geneNames"] += inherited_geneNames
                        existent[taxon]["names"] += [[inherited_taxon, existent[taxon]["names"][-1][1] + inherited_unaCount]]
                        if e_value_enabled:
                            existent[taxon]["eValues"] += inherited_eValues
                            inherited_eValues = []
                        if fasta_enabled:
                            existent[taxon]["fastaHeaders"] += inherited_fastaHeaders
                            inherited_fastaHeaders = []
                        inherited_taxon = ""
                        inherited_geneNames = []
                        inherited_unaCount = 0
                else:
                    if not (taxon in created):
                        created[taxon] = {"taxID": "", 
                                          "children": [],
                                          "unaCount": 0, 
                                          "rawCount": 0, 
                                          "totCount": 0,
                                          "name": name,
                                          "rank": rank, 
                                          "names": [[taxon, -1]], 
                                          "geneNames": []
                        }
                        if fasta_enabled:
                            created[taxon]["eValues"] = []
                        if fasta_enabled:
                            created[taxon]["fastaHeaders"] = []
                    if inherited_unaCount > 0:
                        created[taxon]["unaCount"] += inherited_unaCount
                        created[taxon]["totCount"] += inherited_unaCount
                        created[taxon]["geneNames"] += inherited_geneNames
                        created[taxon]["names"] += [[inherited_taxon, created[taxon]["names"][-1][1] + inherited_unaCount]]
                        if e_value_enabled:
                            created[taxon]["eValues"] += inherited_eValues
                            inherited_eValues = []
                        if fasta_enabled:
                            created[taxon]["fastaHeaders"] += inherited_fastaHeaders
                            inherited_fastaHeaders = []
                        inherited_taxon = ""
                        inherited_geneNames = []
                        inherited_unaCount = 0
            else:
                if j == len(ln) - 1:
                    deleted[taxon] = raw_tax_set[taxon]
                    inherited_taxon = taxon
                    inherited_geneNames = raw_tax_set[taxon]["geneNames"]
                    inherited_unaCount = raw_tax_set[taxon]["rawCount"]
                    if e_value_enabled:
                        inherited_eValues = raw_tax_set[taxon]["eValues"]
                    if fasta_enabled:
                        inherited_fastaHeaders = raw_tax_set[taxon]["fastaHeaders"]
                    lns[i] = lns[i][:j]
                else:
                    lns[i] = lns[i][:j] + lns[i][(j+1):]

    tax_set = {**created, **existent}

    return tax_set, lns

def sort_n_uniquify(lns):
    lns.sort(key=sort_func)

    newlns = copy.deepcopy(lns)
    for i in reversed(range(1, len(lns))):
        if (newlns[i] == newlns[i-1]):
            newlns = newlns[:i] + newlns[i+1:]

    return newlns

def sort_func(ln):
    return ln[-1][1]

def correct_tot_counts(lns, tax_set):
    for i in range(0, len(lns)):
        child = lns[i][-1][1] + " " + lns[i][-1][0]
        tax_set[child]["lnIndex"] = len(lns[i]) - 1
        for j in reversed(range(1, len(lns[i]) - 1)):
            name = lns[i][j][1]
            rank = lns[i][j][0]
            taxon = name + " " + rank
            tax_set[taxon]["totCount"] += tax_set[child]["unaCount"]
            tax_set[taxon]["children"] += [child]
            tax_set[taxon]["lnIndex"] = j
        if child != "root root":
            tax_set["root root"]["totCount"] += tax_set[child]["unaCount"]
            tax_set["root root"]["children"] += [child]
    return tax_set

@app.route('/fetchID', methods = ['POST'])
def fetchID():
    taxName = json.loads(request.data)["taxName"]
    taxid = taxopy.taxid_from_name(taxName, taxdb)
    print("taxid: ", taxid)
    return jsonify({"taxID": taxid[0]})

if __name__ == '__main__':
  app.run(port=5000)
