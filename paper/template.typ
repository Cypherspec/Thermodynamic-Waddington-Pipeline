// Typst template for paper/manuscript.md (pandoc -> typst). Built by paper/build.py.
#let accent = rgb("#0f7d99")
#set document(title: "$title$", author: "$author$")
#set page(
  paper: "us-letter",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(8pt, fill: luma(110))
      Kommisetty · Calibrated thermodynamics of cell fate
      #h(1fr) Preprint
      #v(-4pt)
      #line(length: 100%, stroke: 0.4pt + luma(200))
    ]
  },
)
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, first-line-indent: 0em, spacing: 0.95em)
#show heading.where(level: 1): it => {
  v(0.9em)
  text(13pt, weight: "bold", fill: accent, it.body)
  v(0.25em)
}
#show heading.where(level: 2): it => {
  v(0.55em)
  text(11pt, weight: "bold", it.body)
  v(0.1em)
}
#show figure.where(kind: image): set figure(supplement: [Fig.])
#show figure.caption: it => {
  set text(9pt)
  set par(justify: true)
  align(left, [#text(weight: "bold", fill: accent)[#it.supplement #context it.counter.display(it.numbering).] #it.body])
}
#set figure(gap: 0.7em)
#show figure: it => { v(0.6em); it; v(0.6em) }
#set table(
  stroke: (x, y) => if y == 0 { (top: 0.8pt, bottom: 0.5pt) } else { none },
  fill: (x, y) => if y == 0 { rgb("#eef6f8") },
  inset: (x: 6pt, y: 4.5pt),
)
#show table: set text(8.8pt)
#show table: it => { it; line(length: 100%, stroke: 0.8pt) }
#show link: set text(fill: accent)
#show raw: set text(font: "DejaVu Sans Mono", size: 8.6pt)
#show enum: set text(9pt)
#show enum: set par(spacing: 0.5em)

// Title block
#block(width: 100%)[
  #text(8.5pt, weight: "bold", fill: accent, tracking: 0.08em)[PREPRINT]
  #v(0.3em)
  #text(17pt, weight: "bold")[$title$]
  #v(0.6em)
  #text(11.5pt)[$author$]
  #v(0.15em)
  #text(9pt, fill: luma(90))[$affiliation$ \ Correspondence: $correspondence$ · Code: #link("$code$")]
]
#v(0.6em)
#block(width: 100%, fill: rgb("#f3f8fa"), inset: 11pt, radius: 3pt, stroke: (left: 2.5pt + accent))[
  #text(9.5pt, weight: "bold", fill: accent)[Abstract]
  #v(0.2em)
  #set text(9.6pt)
  $abstract$
  #v(0.3em)
  #text(8.6pt, fill: luma(90))[*Keywords:* $keywords$]
]
#v(0.4em)

$body$
