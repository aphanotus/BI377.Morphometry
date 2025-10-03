# Guidelines for Final Projects

###### tags: `final project`

This document will serve as a guide as you prepare your final projects for BI377. 

## Goal

The purpose of your final project will be to demonstrate your ability to use the methods of morphometry. You are encourage to work in a group of 2 or 3 students.

The project should include:

- identification of an interesting biological question
- obtaining a dataset suitable to address that question, either from publically available sources (including sample datasets provided through the course) or by digitizing images that you acquire
- description of metadata that is relevant to your question, including descriptive statistics, if appropriate
- morphometric analyses to address your question and related issues, with good statistical practice, including reporting of test statistics and p-values
- interpretations and explanation of the results
- ideas for future studies

## How much is enough?

Data gathering is often the limiting step. Publicly available data will be outside of your control, and gathering your own data is an investment of time. Talk to Dave if you plan to digitize your own data. The investment of time will be acknowledged.

Once datasets are in hand, they should be used to address a question that you clearly articulate. However, if your question can be addressed in a straight-forward way, think about supplemental questions that can be addressed using the same data and the tools BI377 has presented. (For example, ANOVA and disparity analysis.) Again, talk with Dave about what scope would be appropriate of your project.

## Products of the BI377 Final Project

You will submit the following products from your project:

- An abstract of no more than 300 words that summarizes the background, question, methods, results and conclusions of your study. This should be suitable for submission to an international conference. 
- Your R code, which should contain annotations identifying the key steps.
- Data files in CSV or TPS format.
- A poster describing the project, in PDF format. 

Posters will be presented at out last class.

## Posters

Posters should be much more detailed than your DataBlitz images. They should convey the background, question, methods, results and conclusions of your study. They should also include citations to the primary literature. (At least 5 references.) Use the citation format of the journal [*Development, Genes and Evolution*](https://www.springer.com/journal/427/submission-guidelines#Instructions%20for%20Authors_Types%20of%20papers). 

[Posters can be printed at the Language Resource Center](https://www.colby.edu/people/offices-directory/language-resource-center/instructional-media-center/poster-project-formatting-and-submission-form/) in Lovejoy hall. Please submit a `pptx` or `pdf` file 48 hours before our session.

### Design

Posters should be 48 x 36 inches, either horizontal or vertical. They include graphics and text and follow the best practices in graphic design and data visualization that our course has discussed. (In general, people usually include too much text. Use text in a supportive role. If something can be shown visually, do that.)


### Title

A concise **Title** should convey the most important conclusions of your study. Titles should state the variable being investigated, what is being measured, and the biological context of the experiment. The best titles will actually convey the conclusions of the study. 

### Author list

Hopefully this is obvious, but please include your names. 

> [!NOTE]
> As an aside, this is sometimes a more profound choice than you may initally think. Once your name appears on a published article in your professional field, that will become the most important search term for people to locate your work. Choose the inclusion of middle names carefully. For people with common names, there is a danger of confusion. Also, if your name comes from a culture where the [surname comes first](https://en.wikipedia.org/wiki/Surname#Order_of_names), will that be obvious to readers? A trend in the sciences is to couple names with a unique identiying number such as [ORCID](https://orcid.org/). (e.g. [0000-0002-2776-2158](https://orcid.org/0000-0002-2776-2158)).

### Background

Do not include your abstract on your poster. 

However, a poster should include key background information necessary for readers to understand the purpose and interpretation of your experiment. You want to describe the issue or question being addressed by your study. You will need to provide background from the [published literature](https://hackmd.io/@aphanotus/BI376_writing_guide#Finding-sources-from-the-published-literature) for this assignment. 

The Introduction should be organized so that it moves from general information to more specific information. It should begin by emphasizing the importance of the topic you're studying. Details related to your specific study system and experiments should come later in the section.

#### Hypothesis

At some point in your Introduction, you should make clear what hypothesis your study will test. Ideally, this is done early, such as the end of the first paragraph. Remember, a hypothesis is a general statement about the process you predict to explain a previously observed pattern. After stating your hypothesis, you should describe how it will be tested.

> [!TIP]
> * A common pitfall is that students present predictions as though they were a hypothesis. But a hypothesis is meant to be a general statement about the biological process being studied. Your predictions are the results you expect from the experiment, if your hypothesis is actually true. 
> * When you state your hypothesis, try to avoid using "hypothesize" as a verb. It's not wrong, but it is inelegant. Some readers (including your instructor) find it grating!

It is not necessary for your results to support your hypothesis. It is possible for you to produce an excellent project, in which you present a hypothesis whose predictions are not met by your results. Scientific articles get a lot of attention if they challenge or overturn conventional thinking in their field. So authors will often pose their hypothesis based on that conventionial thinking (and supported by the background they present in the Intro). If the results don't match the predictions, the Discussion can then speculate on how and why everyone's previous understanding was incomplete. In some cases people even pose a hypothesis that's so weak, a skeptical reader might call it a "strawman" -- one that's easily knocked down. It's okay to do this in our class if you need to!

### Materials & Methods

This section should be brief and describe each phase of your study from data gathering through the analysis.  If your data are based on published studies, be sure to include references to those sources. Be sure to include a detailed description of anatomical landmarks used in GMM. 

Any R packages used in the analysis should be named in the Methods section too. Packages should also be cited. Many will be described in published papers; some will just cite a web site. If you look at the help page for a function in a package, it should tell you how to cite it. Otherwise, you can use the R function `citation`, as below.

```{r}
citation("geomorph")
citation("borealis")
citation("Momocs")
```

### Results & Discussion

For this assignment you are welcome to combine results and their interpretation into one section, but you may also choose to keep them separated. Sub-section headings should communicate results in a short phrase. The results should be clearly organized and summarized in writing for a reader to judge. Provide the details that are necessary to reach each conclusion. One way to do this is to state the expected results. If you are combining Results & Discussion, talk through the logic that brings a reader to your conclusion. 

> [!TIP]
> * Reference your figures parenthetically, like this (Fig. 1). 
> * In a reference to a figure or source, never include the word “See”. -- Don't do this (See Fig. 1).
> * Don't ignore results just because they don't match what you expect or want to see. 

For each analytical or statistic test, you should explain its purpose and the data included in the test. Provide enough of an explanation to allow the readers to understand the test and interpret its results.

Keep in mind that most details should be covered in the [Materials & Methods](https://hackmd.io/@aphanotus/BI376_writing_guide#Materials-amp-Methods) section. Remember that to properly [report the results of a statistcal test](https://hackmd.io/@aphanotus/bi376f20_LabExercise2#Reporting-the-results-of-statistical-tests) you must name the test, include the test statistic and the p-value.

In each sub-section (or in a distinct Discussion section) provide your interpretation of the results of the experiment. What do the results mean? 

#### Quantitative data

Typically it is best to also report quantitative data in some sort of graph. Whenever you want to say two things are different, that assertion should be supported by an appropriate [statistical test](https://hackmd.io/@dts8RULgQqi0n0PPDKh7JQ/rJmbl8Q2S), unless there is [complete separation](https://en.wikipedia.org/wiki/Separation_(statistics)) (no overlap) in the data. In addition to the exercises covering the statistics of shape analysis, a [guide to statistics](https://hackmd.io/@dts8RULgQqi0n0PPDKh7JQ/rJmbl8Q2S) also includes detailed instructions for performing basic statistical tests in R. You can also talk to Dave!

> [!WARNING]
>
> * "Data" are plural. So you should write, "Our data show..." rather than "Our data shows..." (The singular form of "data" is "datum.")

#### Figures

Figures should be used to present images or graphical data. Think carefully about how best to organize the visual impact of image data. Label the important parts of anatomical images. Think about graphic design for the best communication of your data. For example, crop blank space from the periphery of pictures of animals. Images that need to be compared to one another should be included as separate panels within the same figure. 

Many students find PowerPoint a convenient program in which to compose figures. [ImageJ](https://imagej.net/Welcome) is free software that can be used to crop or rotate images, adjust brightness and contrast, along with many other image manipulations.

Always include a **legend** for figures and tables that gives a brief sentence summarizing the figure's meaning—not its content. Think of this as the title to the figure. ("Figure 1. *engrailed* is expressed in the posterior compartment of the *Drosophila* wing.") Individual panels within the figure should be lettered, and important content should be labeled. Since space may be tight, abbreviations are okay, but define them in the figure legend. After the figure title, explain the content of each panel. 

> [!TIP]
> * Keep multiple specimens oriented in the same direction!
> * Don’t make figures or tables the actors in your writing. (For example, avoid statements like "Figure 1 shows that our PCR was successful." or "Our results are shown in Table 1.") Just cite results parenthetically. (e.g. "Misexpression of human *FoxP2* in the fly brain resulted in increased language ability (Fig. 2).") 

#### Tables

In rare cases, a small table may be okay on a poster, but in general they are not well suited to this communication format and should be avoided on posters. Instead, think of ways to communicate the information graphically. (In a poster a bulleted list is okay too!) If you just use a table, be sure to give headings to each column and include units where appropriate. 

### Conclusions

Use the last section of your poster to return to your hypothesis and state the overall conclusions unambiguously. You should also put your study in a broader context and point out the significance of your work. If your experiment was unable to adequately address the question, you can offer an explanation of why that might be so. Discuss what might be done in the future to better accomplish your goal or to further investigate the topic.

> [!TIP]
> * Don't undersell your results!
> * Don't  present a litany of problems. 
> * Don't ascribe unexpected results to human error. If needed be specific about  unexpected problems that prevented you from getting a clear result.

### Acknowledgements 

This section is optional on a poster. In the Acknowledgements you can give credit to anyone who may have helped you with experiments or analysis. It is also appropriate to thank anyone who provided significant help with your writing process.

### Literature Cited

A section of references should list the full citations of any source referenced in the text. Follow the [citation format guidelines](https://www.nlm.nih.gov/bsd/uniform_requirements.html) recommended by the International Committee of Medical Journal Editors and the US National Institutes of Health.

#### Finding sources from the published literature

Search the literature for articles that provide background information on your topic. In a busy field this can be daunting, but it is essential. The good news is that it's not necessary to read everything! (At least while you're in college.) For any topic you will probably find hundreds of articles. Choose your search terms carefully to provide a focused, manageable number of titles. Use titles and abstracts to filter down to the articles that are most useful to you. For example, scan the titles of about 25 articles. Choose the most relevant and carefully read the abstracts of about 10. From those, read only a handful of articles in detail. If what your project touches on multiple topics, you should repeat this process for each search-able set of terms you have in mind. As you read, be open to searching more if the articles you read suggest there are other areas you should think about. 

> [!TIP]
>
> One important way to find useful articles is to search "backwards" and "forwards". Go backward in time by looking at the references within articles to find other useful ones. Then search forward in time by using citation databases to find [later articles that cite](https://www.ncbi.nlm.nih.gov/pubmed?linkname=pubmed_pubmed_citedin&from_uid=15483612) those you've already found. A professor or the [science librarians](http://www.colby.edu/olin/) can help you use these databases.

Many search tools exist for the literature. [Google Scholar](https://scholar.google.com/) will provide fast results, but search fields are not very customizable and links to journals may require payment for access. Colby subscribes to several databases that allow you to search for keywords, journals, authors, years, etc. The links below will provide you access to all journals to which the college has active subscriptions.

**Literature Search Tools** 

- [PubMed](https://colby.idm.oclc.org/login?url=http://www.ncbi.nlm.nih.gov/entrez/query.fcgi) - The US public database of biomedical relevant literature. 
- [SCOPUS](https://colby.idm.oclc.org/login?url=https://www.scopus.com/) - A commercial search engine for academic and scholarly literature.
- [JStor](https://colby.idm.oclc.org/login?url=https://www.jstor.org) - A database covering less medically-oriented sciences, including ecology and organismal biology, as well as the humanities and social sciences.

> [!TIP]
> [ILL](https://colby-illiad-oclc-org.colby.idm.oclc.org/illiad/illiad.dll) - If you have difficulty getting access to an article, you can request an "inter-library loan". Complete the form, and library staff at Colby will find the article or book anywhere in the world and provide you with a PDF (for articles) or a hard copy (for books). 

As you search, keep notes for yourself on the ways you might use these articles in your paper. While the words "annotated bibliography" might make you grown, thinking back to high schools English, in our class, these notes are not graded. This is just good practice, and it can help you stay organized and ultimately save you time.

> [!TIP]
> Your notes might look like this. You can also paste URLs into your notes to link back to the articles.
>
> * Carroll 1995 *Nature*. [https://www.ncbi.nlm.nih.gov/pubmed/7637779](https://www.ncbi.nlm.nih.gov/pubmed/7637779) -- This is a review that explains how homeotic genes are conserved across animals. Clusters of them may evolve through gene duplication. Good ref for these background facts.
> * Halder, Callaerts & Gehring 1995 *Science* [https://www.ncbi.nlm.nih.gov/pubmed/7892602](https://www.ncbi.nlm.nih.gov/pubmed/7892602) -- Using Gal4/UAS they ecotopically expressed *eyeless/Pax6* in the legs of flies and got eyes to develop on the leg! A good ref for the  importance of *Pax6* in eye development: necessary and sufficient. Might also serve as a ref showing the power of Gal4/UAS misexpression.
> * Yamamoto, et al. 2004 *Nature* [https://www.ncbi.nlm.nih.gov/pubmed/15483612](https://www.ncbi.nlm.nih.gov/pubmed/15483612) -- Cite this one to show Hedgehog is important for eye development in vertebrates. Cool experiment showing that applying a drug to inhibit Hh partly rescues development of eyes in "eyeless" cave fish.



## Advice

### Academic Integrity

The scientific community and the Colby community value [intellectual honesty](http://www.colby.edu/academicintegrity/the-colby-affirmation/). The work you present must be your own, and you should always give credit when information or ideas were developed by others. Plagiarism includes misrepresenting any part of a written document as your own if it isn't. Uncredited text of any length (cutting and pasting as well as paraphrasing) also violates academic integrity. You may be tempted to think that a source you read says something better than you ever could. But you are writing *your* ideas about their work. So, describe it in *your* words and give them credit for the idea or information with a citation. Colby's libraries have more resources for [avoiding plagiarism](http://libguides.colby.edu/avoidingplagiarism).

The only circumstances where citations aren't necessary, are for "common knowledge" in your field. For example, it's not necessary to cite [Darwin (1859)](https://www.gutenberg.org/files/1228/1228-h/1228-h.htm) for the theory of natural selection. However, if you're uncertain, it's always best to be cautious and err on the side of making unnecessary citations rather than being called out for plagiarism. 

> [!CAUTION]
>
> * Never include direct quotes from your sources, even if they *are* cited. Why? This practice reflects the philosophy that we should value the information, independent of the personality of its authors. We credit them, but directly quoting is seen as too reverent.


### Proofreading

After you have a draft, proofread your paper—twice! First, read it over to look at the narrative flow. Do ideas flow logically from one point to the next? Do you have sensible headings? Someone else may be better for this purpose. You may be emotionally attached to your own words. So ask a friend to be merciless, and take her advice when she says to cut that sentence. Next, edit for grammar, spelling, and the formatting of references, chemical symbols, or mathematically notation. (Don't get bogged down in these details early, just to have them cut or overhauled later because the sentence makes no sense.) Give yourself the time to write. Life is always busy. But when you use time to work on writing, it will improve. 

### Voice and tense

In the Introduction you'll present background information that is mostly accepted as true. These facts should be conveyed in the present tense. ("Frogs accomplish gas exchange via diffusion across the skin.") The results of your experiment should be written in the past tense. Avoid the [passive voice](https://webapps.towson.edu/ows/activepass.htm). ("Cells were cultured overnight." -- This point of grammar has seen a lot of [debate](http://advice.writing.utoronto.ca/types-of-writing/active-voice-in-science/) in the last several decades, but most people now favor the use of a clearer, active voice.)

Scientific writing does not emphasize the author. But you will probably still use personal pronouns like "I" or "we", but be sparing. ("We have conduced a detailed study of acorns, and our results reveal...")

> [!TIP]
> * Don't use the word "gender" in reference to animals or plants. Gender is a human idea that's irrelevant to other organisms. Biologically speaking, your refering to "sex" when you refer to an individual animal's status as male or female (or [intersex](https://doi.org/10.1016/j.ydbio.2011.09.026)!). 

### Species names

Species names always appear in italics, such as the fruit fly *Drosophila melanogaster*. The first mention of a species in an article should use the full name (genus and specific epithet). Properly, subsequent mention of the same species should abbreviate the genus and use the specific epithet, as in *D. melanogaster*. While it is not strictly correct, molecular biologists often simply use the genus name for subsequent mentions of the same organism, for example referring again to *Drosophila*. Authors may also interchange the scientific and common names (flies and *Drosophila*). 

Higher taxonomic groups, like family and phylum are treated as proper names and capitalized. However, the names can often be used as non-italics non-Latin names. ("Diptera is the order containing flies." vs. "Insecta is the class containing dipterans.") Taxonomists will insist that only "natural groups" (monophyletic groups) be given Latin names. Since it is now accepted that protists are not monophyletic, "Protista" should no longer be used. Some group names can also be used as adjectives. ("Heteroptera includes the cicadas, seed bugs and bed bugs. These insects share a specialized heteropteran mouthpart anatomy.")

For more info and stories related to scientific names check out the blog [*Scientist Sees Squirrel*](https://scientistseessquirrel.wordpress.com/category/latin-names/).

### Abbreviations

Overuse of abbreviations and uncommon jargon can annoy readers. If something will be referred to many times and it has an abbreviation or technical (jargon) term, you are entitled to use it, but define it at your first use. ("The occurrence of adults in some insect species with long and short wings is known as brachyptery." "The enzyme $\beta$-galactosidase ($\beta$-gal) can break-down the organic substrate molecule 5-bromo-4-chloro-indolyl-$\beta$-D-galactopyranoside (X-Gal).") Some abbreviations are so common that no definition is needed, like DNA or PCR. If you're uncertain, ask yourself if friends in the biology major would be sure to know the abbreviation.

There are some common abbreviations and foreign phrases that appear in science that deserve to be defined (Table 1).

**Table 1: Common abbreviations and non-English phrases used in biology.**

| **Notation**             | **Definition**                                               |
| ------------------------ | ------------------------------------------------------------ |
| e.g.                     | "for example" (Latin); don't use with "etc." since it's redundant    |
| etc.                     | "and the rest" (Latin)                                               |
| et al.                   | "and others" (Latin); usually meaning "and coworkers" in an author list |
| *in vivo*                  | "in life" (Latin); in a live organism                                |
| *in vitro*                | "in glass" (Latin); in an artificial system                          |
| *in silico*                | "in a computer simulation" (Latin)                                   |
| i.e.                     | "that is" (Latin); used to start a complete list; incomplete lists start with "e.g." |
| *senu stricto*             | "in the strictest sense" (Latin); the narrowest definition of the subject |
| *senu lato*                | "in the loose sense" (Latin); the broadest definition of the subject |
| locus (plural: loci)     | In genetics, a locus (Latin for location) refers to a position on a chromosome. It can refer to a gene, or any other large genetic element, like a transposon. People tend not to use the term when you get down to the level of nucleotide sequences. Don't confuse the singular (locus) and plural (loci -- pronounced "low sigh") |
| milieu                   | A French word meaning "surroundings". Often used by chemists to mean other compounds in solution, or by cell biologists to refer to other things in a cell, which are not immediately under consideration. |
| anlage (plural: anlagen) | From German, most literally meaning "primordium". However, it can be used to refer to any "developmental field" or "territory" at any stage of development. |

These guidelines include the basics, and some of the relatively easy things you can do to give your writing more polish. The best thing you can do to develop your skills in science writing is to practice and keep reading the primary scientific literature.

------



## Checklist & Rubric

#### Higher-order concerns

- [ ] Adopt a tone and style appropriate to an audience of peer scientists
- [ ] Write concisely, avoiding extraneous words or phrases
- [ ] Write with scientific accuracy, especially for terminology and essential biological concepts 
  + Write to emphasize the biological entities rather than figures or citations
- [ ] Structure the poster in a logical progression 
- [ ] Demonstrate an understanding of the scientific method
  + Distinguish among hypotheses and predictions
  + Qualify speculations and interpretations as such
  + Support assertions with evidence
- [ ] Present results precisely and accurately
  + Explain the purpose or significance of experiments
  + Effectively convey information visually (e.g. choice of figure type, clear labels)
- [ ] Analyze results 
  + Summarize patterns in the data that support or refute a hypothesis
  + Use of statistical tools, as appropriate, to identify patterns or differences among groups
- [ ] Discuss findings in a way that considers alternative explanations (i.e. avoid confirmation bias)
- [ ] Synthesize results with prior knowledge

#### Format

- [ ] Title is appropriate
- [ ] Author list is included
- [ ] Abstract summarizes key points from all sections
- [ ] Background explains the importance and purpose of the study
- [ ] Materials & Methods section succinctly describes key components of the analyses
- [ ] Materials & Methods cite the sources of key samples, previously published data and analytical methods
- [ ] Analytical code is included
- [ ] Results & Discussion section mentions procedures, but not details of methods
- [ ] Statistical tests are named and presented with values of the test statistic, *p*-value, and interpreted in a biological context
- [ ] Figures are well-designed to communicate the data
- [ ] Figures are numbered in the order in which they're discussed in the text
- [ ] Figures and tables have descriptive legends
- [ ] Discussion interprets the results
- [ ] Discussion / Conclusions connects to the background information to highlight broader significance
- [ ] Citations are given wherever appropriate, throughout the paper
- [ ] Reference list gives full information for all citations
- [ ] Reference list is properly formatted

#### Graphic design style

- [ ] Shows the data without distortion 
- [ ] A visual strategy organizes data
- [ ] Encourages thought about the subject (not the method or presentation)
- [ ] Design organizes the flow of the study in a logical direction
- [ ] Design avoids unnecessary distractions
- [ ] Design is eye-catching

#### Writing style

- [ ] Species names are italics
- [ ] Terminology is used properly
- [ ] Uncommon abbreviations are define before use
- [ ] Measurable values have units

#### Figures

- [ ] Descriptive legends
- [ ] Panels are labeled
- [ ] Consistent, appropriate framing
- [ ] Consistent, appropriate orientation
- [ ] Important results are clearly evident

#### File submissions

BI377 final projects should include:

- [ ] Poster in `pdf` format
- [ ] Landmark data (if your project involved new data gathering) in `tps` format
- [ ] Metadata in `csv` format
- [ ] Analytical code in `R` script or `Rmd` (R markdown) format




---

**Quick Links** | [BI377 Moodle](https://moodle.colby.edu/course/view.php?id=33250) | [BI377 HackMD](https://hackmd.io/@ColbyBI377/landingpage) | [rstudo.colby.edu](https://rstudio.colby.edu/)

---
