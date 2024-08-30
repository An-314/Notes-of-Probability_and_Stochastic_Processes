#import "@local/mytemplate:1.0.0": *

#show: project.with(
  title: "概率论与随机过程",
  template: "book",
  authors: (
    "Anzreww",
  ),
  time: "甲辰春夏于清华园",
  contents: true,
  content_depth: 3,
)

#pagebreak(weak: true)

#include "chap1.typ"

#pagebreak(weak: true)

#include "chap2.typ"
