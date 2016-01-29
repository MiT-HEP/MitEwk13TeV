import copy
import subprocess

class LatexDocument:
    def __init__(self):
        self.texFileName = "results.tex"
        self.figDir = "plots/"
        self.texFile = None

    def openDocument(self):
        self.texFile = open(self.texFileName, "w")
        print >>self.texFile, "\\documentclass[a4paper,12pt]{article}"
        print >>self.texFile, "\\usepackage{subfigure}"
        print >>self.texFile, "\\usepackage{amssymb}"
        print >>self.texFile, "\\usepackage{multirow}"
        print >>self.texFile, "\\usepackage{rotating}"
        print >>self.texFile, "\\usepackage{url}"
        print >>self.texFile, "\\renewcommand{\\topfraction}{1.0}"
        print >>self.texFile, "\\renewcommand{\\bottomfraction}{1.0}"
        print >>self.texFile, "\\renewcommand{\\textfraction}{0.0}"
        print >>self.texFile, "\\newlength{\\dinwidth}"
        print >>self.texFile, "\\newlength{\\dinmargin}"
        print >>self.texFile, "\\setlength{\\dinwidth}{21.0cm}"
        print >>self.texFile, "\\textheight23.5cm \\textwidth16.0cm"
        print >>self.texFile, "\\setlength{\\dinmargin}{\\dinwidth}"
        print >>self.texFile, "\\setlength{\\unitlength}{1mm}"
        print >>self.texFile, "\\addtolength{\\dinmargin}{-\\textwidth}"
        print >>self.texFile, "\\setlength{\\dinmargin}{0.5\\dinmargin}"
        print >>self.texFile, "\\oddsidemargin -1.0in"
        print >>self.texFile, "\\addtolength{\\oddsidemargin}{\\dinmargin}"
        print >>self.texFile, "\\setlength{\\evensidemargin}{\\oddsidemargin}"
        print >>self.texFile, "\\setlength{\\marginparwidth}{0.9\\dinmargin}"
        print >>self.texFile, "\\marginparsep 8pt \\marginparpush 5pt"
        print >>self.texFile, "\\topmargin -42pt"
        print >>self.texFile, "\\headheight 12pt"
        print >>self.texFile, "\\headsep 30pt \\footskip 24pt"
        print >>self.texFile, "\\parskip 3mm plus 2mm minus 2mm"
        print >>self.texFile, "\\begin{document}"

    def addSection(self, section):
        print >>self.texFile, "\\section{"+section+"}";

    def insertFigure(self, name, caption):
        print >>self.texFile, "\\begin{figure}[!htbp]";
        print >>self.texFile, "\\begin{center}";
        print >>self.texFile, "\\begin{tabular}{lr}";
        print >>self.texFile, "\\includegraphics[width=14cm]{"+self.figDir+"/"+name+"}";
        print >>self.texFile, "\\end{tabular}";
        print >>self.texFile, "\\end{center}";
        print >>self.texFile, "\\caption{"+caption+"}";
        print >>self.texFile, "\\end{figure}";
        print >>self.texFile, "\\newpage";

    def closeDocument(self):
        print >>self.texFile, "\\end{document}";
        self.texFile.close()

    def compile(self):
        fileNameCopy = copy.copy(self.texFileName)
        subprocess.call(["latex", self.texFileName])
        subprocess.call(["dvipdf", fileNameCopy.replace(".tex", ".dvi")])
