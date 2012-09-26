import datetime
import time
tz = time.tzname[0]
now = datetime.datetime.now()

intro = r"""
<html>

<head>
<title>Arb documentation</title>

<style>
body { font-family: sans-serif; font-size: 90%; }
dt { color:#348; font-weight: bold; }
dd p { margin-top: 0.4em; }
</style>

<script type='text/x-mathjax-config'> 
MathJax.Hub.Config({
    extensions: ["tex2jax.js"],
    jax: ["input/TeX", "output/HTML-CSS"],
    MMLorHTML: { prefer: "HTML" },
    tex2jax: {
      inlineMath: [ ['$','$'], ["\\(","\\)"] ],
      displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
      processEscapes: true
    },
    "HTML-CSS": { availableFonts: ["TeX"] }
  });</script> 
<script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script> 

<script type='text/javascript'> 
  newContainer = document.createElement('span');
  newContainer.style.setProperty("display","none","");
  newNode = document.createElement('script');
  newNode.type = "math/tex";
  newNode.innerHTML = '\\newcommand{\\NN}{\\mathbb{N}}\n\\newcommand{\\ZZ}{\\mathbb{Z}}\n\\newcommand{\\QQ}{\\mathbb{Q}}\n\\newcommand{\\FF}{\\mathbb{F}}\n\\newcommand{\\KK}{\\mathbb{K}}\n\\newcommand{\\RR}{\\mathbb{R}}\n\\newcommand{\\CC}{\\mathbb{C}}\n\\newcommand{\\abs}[1]{\\left|#1\\right|}';
  newContainer.appendChild(newNode);
  document.body.insertBefore(newContainer,document.body.firstChild);
</script> 

</head>

<body>

<h1>Arb documentation</h1>

<p><i>Last updated: %NOW%</i></p>

<h2>Contents</h2>

"""

intro = intro.replace("%NOW%", "%i-%02i-%02i %02i:%02i:%02i %s" % (now.year, now.month, now.day, now.hour, now.minute, now.second, tz))

outro = r"""

<p>Wishes or bug reports? Send any comments to <a href="mailto:fredrik.johansson@gmail.com">fredrik.johansson@gmail.com</a></p>

</body>
</html>
"""

def section_key(title, section):
    def stringify(s):
        t = ""
        for c in s:
            if c.isalpha():
                t += c
            else:
                t += "-"
        return t
    return stringify(title) + "-" + stringify(section)

def parse(string, markup):
    if not markup:
        return string
    lines = string.splitlines()
    # Collect paragraphs
    chunks = []
    chunk = []
    for line in lines:
        if not line.strip():
            if chunk:
                chunks.append(chunk)
                chunk = []
        else:
            chunk.append(line)
    doc = []
    heading = ""
    description = []
    in_heading = False
    in_description = False
    chunks.append(["***"])
    for chunk in chunks:
        if len(chunk) == 1 and chunk[0].startswith("***"):
            if in_heading:
                in_heading = False
                doc.append(("H", heading.lstrip()))
                heading = ""
            else:
                in_heading = True
                if in_description:
                    doc.append(("DESCR", description))
                    description = []
                    in_description = False
        else:
            if chunk[0].startswith(" "):
                if in_heading:
                    heading += " ".join(chunk)
                else:
                    in_description = True
                    description.append(" ".join(chunk))
            else:
                if in_description:
                    doc.append(("DESCR", description))
                    description = []
                    in_description = False
                doc.append(("DEF", " ".join(chunk)))
    return doc

def write(docs, file):
    file.write(intro)
    file.write("<ul>\n");
    # Write TOC
    for doc, title in docs:
        file.write("<li>\n");
        file.write('<a href="#%s">%s</a>\n' % (section_key(title, ""), title));
        file.write("<ul>\n");
        if isinstance(doc, list):
            for entry in doc:
                if entry[0] == "H":
                    file.write('<li><a href="#%s">%s</a></li>\n' % (section_key(title, entry[1]), entry[1]))
        file.write("</ul></li>\n");
    file.write("</ul>\n");
    file.write("<hr />\n")
    # Write each section
    for doc, title in docs:
        in_list = False
        file.write('<h2><a name="%s">%s</a></h2>\n' % (section_key(title, ""), title));
        if isinstance(doc, list):
            for entry in doc:
                if entry[0] == "H":
                    if in_list:
                        in_list = False
                        file.write("</dl>\n")
                    file.write('<a name="%s"><h3>%s</h3></a>\n' % (section_key(title, entry[1]), entry[1]))
                if entry[0] == "DEF":
                    if not in_list:
                        in_list = True
                        file.write("<dl>\n")
                    #file.write("<dt><tt>%s</tt></dt>\n" % entry[1])
                    file.write("<dt>%s</dt>\n" % entry[1])
                if entry[0] == "DESCR":
                    if in_list:
                        file.write("<dd>\n")
                    for para in entry[1]:
                        file.write("<p>%s</p>\n" % para)
                    if in_list:
                        file.write("</dd>\n")
            if in_list:
                file.write("</dl>\n")
        else:
            file.write(doc)
        file.write("<hr />\n")
    file.write(outro)

docs = [
  ("introduction.txt", "Introduction", False),
  ("setup.txt", "Setup", False),
  ("fmpr.txt", "fmpr.h (floating-point arithmetic)", True),
  ("fmprb.txt", "fmprb.h (real ball arithmetic)", True),
  ("fmprb_poly.txt", "fmprb_poly.h (polynomials of real balls)", True),
  ("fmprb_mat.txt", "fmprb_mat.h (matrices of real balls)", True),
  ("history.txt", "History", False),
  ("credits.txt", "Credits", False),
]

write([(parse(open(doc).read(), markup), title) for \
    (doc, title, markup) in docs], open("doc.html", "w"))

