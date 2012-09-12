import datetime
now = datetime.datetime.now()

intro = r"""
<html>

<head>
<title>Arb</title>

<style>
body { font-family: sans-serif; font-size: 90%; }
dt { color:#348; }
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

<!--
<script type='text/javascript'> 
  newContainer = document.createElement('span');
  newContainer.style.setProperty("display","none","");
  newNode = document.createElement('script');
  newNode.type = "math/tex";
  newNode.innerHTML = '\\newcommand{\\NN}{\\mathbb{N}}\n\\newcommand{\\ZZ}{\\mathbb{Z}}\n\\newcommand{\\QQ}{\\mathbb{Q}}\n\\newcommand{\\FF}{\\mathbb{F}}\n\\newcommand{\\KK}{\\mathbb{K}}\n\\newcommand{\\RR}{\\mathbb{R}}\n\\newcommand{\\CC}{\\mathbb{C}}\n\\newcommand{\\abs}[1]{\\left|#1\\right|}';
  newContainer.appendChild(newNode);
  document.body.insertBefore(newContainer,document.body.firstChild);
</script> 
-->

</head>

<body>

<h1>Arb</h1>

<p><i>Last updated: %NOW%</i></p>

<p>Arb is an experimental C library implementing arbitrary-precision floating-point ball arithmetic, written by Fredrik Johansson &lt;<a href="mailto:fredrik.johansson@gmail.com">fredrik.johansson@gmail.com</a>>. The git repository is <a href="https://github.com/fredrik-johansson/arb/">https://github.com/fredrik-johansson/arb/</a></p>

<p>Ball arithmetic, also known as mid-rad interval arithmetic, is an extension of floating-point arithmetic in which an error bound is attached to each variable. This allows doing rigorous computations over the real numbers, while avoiding the overhead of traditional (inf-sup) interval arithmetic at high precision, and eliminating much of the need for time-consuming and bug-prone manual error analysis.</p>

<p>At the moment, Arb contains:</p>

<ul>
<li>A module (fmpr) for correctly rounded arbitrary-precision floating-point arithmetic. Arb numbers have a few special features, such as arbitrary-size exponents (useful for combinatorics and asymptotics) and dynamic allocation (facilitating implementation of hybrid integer/floating-point and mixed-precision algorithms).</li>
<li>A module (fmprb) for real ball arithmetic, where a ball is implemented as a pair of fmpr numbers.</li>
<li>Functions for fast high-precision computation of some mathematical constants, based on ball arithmetic.</li>
<li>A rudimentary module (fmprb_poly) for polynomials or power series over the real numbers, implemented using balls as coefficients, with fast polynomial multiplication.</li>
</ul>

<p>Planned features include: transcendental functions, more extensive polynomial functionality, matrices, and complex balls (and polynomials and matrices thereof).</p>

<p>Arb uses <a href="http://mpir.org/">MPIR</a> and <a href="http://flintlib.org/">FLINT</a> for the underlying integer arithmetic. The code conventions borrow from FLINT, and the project might get merged back into FLINT when the code stabilizes in the future. It also uses <a href="http://mpfr.org">MPFR</a> for some fallback code and for testing purposes. The current version of Arb implements most of its floating-point arithmetic naively using high-level FLINT types. The speed at low precision is far from optimal, and the memory management can sometimes be wasteful. The internals will be rewritten in the future to fix the inefficiencies, which eventually should make Arb ball arithmetic about as fast as mpz or mpfr arithmetic at any precision.</p>

<p><b>Warning</b>: as this is an early version, any part of the interface is subject to change! Also be aware that there are known and unknown bugs.</p>


<h2>Contents</h2>

"""

intro = intro.replace("%NOW%", "%i-%02i-%02i %02i:%02i:%02i" % (now.year, now.month, now.day, now.hour, now.minute, now.second))

outro = r"""


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

def parse(string):
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
        for entry in doc:
            if entry[0] == "H":
                file.write('<li><a href="#%s">%s</a></li>\n' % (section_key(title, entry[1]), entry[1]))
        file.write("</ul></li>\n");
    file.write("</ul>\n");
    # Write each section
    in_list = False
    for doc, title in docs:
        file.write('<h2><a name="%s">%s</a></h2>\n' % (section_key(title, ""), title));
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
    file.write(outro)

docs = [
  ("fmpr.txt", "fmpr (floating-point arithmetic)"),
  ("fmprb.txt", "fmprb (real ball arithmetic)"),
  ("fmprb_poly.txt", "fmprb_poly (polynomials of real balls)"),
]

write([(parse(open(doc).read()), title) for (doc, title) in docs],
    open("doc.html", "w"))

