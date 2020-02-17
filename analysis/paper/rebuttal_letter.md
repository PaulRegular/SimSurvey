Daniel E. Duplisea, PhD  
Fisheries and Oceans Canada  
850 Route de la Mer  
Mont‐Joli, Quebec  
G5H 3Z4, Canada

2020-02-17

Dear Dr. Duplisea,

Thank you for considering another revision of manuscript
PONE-D-19-26904R1, **“`SimSurvey`: an `R` package for comparing the
design and analysis of fisheries surveys by simulating
spatially-correlated fish stocks”** by Paul M. Regular, Gregory J.
Robertson, Keith P. Lewis, Jonathan Babyn, Brian Healey and Fran
Mowbray. We are also greatful for your detailed suggestions and we have
made every effort to do justice to the changes you reccomend. Most
importantly, we hope we have added sufficient content to the core of the
paper elevate the manuscript from soley a software manual to a primary
scientific publication. Though the how-to approach remains, we now see
that describing some of the case study results in the core of the
manuscript makes it more interesting and it adds another tangible reason
for prospective users to learn how to use the package. Please see below
for more details on the changes we made in response to your suggestions.

We submit this revised manuscript for your consideration and look
forward to your decision.

Sincerely,

Paul Regular  
Fisheries and Oceans Canada  
Northwest Atlantic Fisheries Center  
80 East White Hills, St. John’s, NL  
A1C 5X1, Canada  
E-mail:
<a href="mailto:Paul.Regular@dfo-mpo.gc.ca" class="email">Paul.Regular@dfo-mpo.gc.ca</a>  
Phone: (709) 772-2067

------------------------------------------------------------------------

Additional Editor Comments (if provided):

Dear Paul, you have made lots of excellent changes that I think help
with the use of the package and responded to the many specific comments
of reviewers which has no doubt corrected many technical issues and
clarifications. The one thing I would say that you have not really done
is get at the deeper research merits of this work beyond introducing a
new piece of software. I think paper, as it stands, lacks a larger
context and content which is important for the primary publication. My
diagnosis for why this is is that the manuscript does not conform very
well to more typical scientific reports (Intro, M&M, Result, Discussion)
which can make it difficult for readers to find the larger scientific
merits of the work. It is useful of course for those who already
understand the merits of this kind of work but this work is for primary
publication and it needs to appeal more to the former than the latter
group. There is a very “how-to” feel to it (e.g. line 64 “In this
section”) which I think detracts from getting at the larger purpose of
the work.

I would really like you to address this issue of moving it from a
software manual to a primary scientific publication. I do not think it
should involve that much work but there will be some restructuring of
sections as well as places to put in content and bring out conclusions.
Here are my suggestions for this:

Try to follow a more traditional paper structure. This will help readers
and it likely will also make it clearer for you on how you can inject
content into the paper to move it beyond the software manual approach:

Introduction:

-   You need to talk about wider issues and examples. e.g. examples of
    when poor survey design meant that scientific questions could not be
    properly addressed when good survey design meant they were. Examples
    of when good survey design allowed researchers to address needs that
    were unanticipated at the time it was designed. i.e. you need to
    build a better case (not just cost) of why survey design is really
    important and it is most powerful to do this with examples. You
    should try to bring in ideas related to ecosystem and climate
    changes and being able to track communities. Perhaps bring in
    species at risk ideas and tracking decline, you could bring in ideas
    related to MSC certification. These are just examples of specifics
    but you get the idea: the Introduction needs to have more general
    information outlining in both a broad and specific sense why we
    should be concerned about this.

-   Keep in mind an educated reader who may not be in fisheries but is
    interested in why anyone should care about this or could, for
    example, be interested in surveying say caribou or songbirds or
    something outside of fisheries but where many of the motivation and
    concepts may be similar and they are doing a more general literature
    search before designing their own survey.

Methods:

-   This will start at your “Model Structure” section. You should put a
    higher level heading just before that called Methods.

-   I can see that it is hard to separate your Methods from Results.
    Your results are really in the section “Using SimSurvey” I would not
    be opposed to putting this as a separate Results section but I leave
    it to you to decide. Essentially what you have is a case study as an
    example of how it works and therefore specific results are less
    interesting than how you got there with the package. You might title
    it something like “Results: running a SimSurvey simulation”. This
    section also has a lot of content without a lot of explanation. for
    example, you have several large and complicated figures in a row.
    You should discuss not only what figures are meant to offer in the
    package but also what they mean. So for example on line 473 you
    refer to three figures (5,6,7) but you offer little interpretation
    of those figures just why you can make them which is another example
    of the limitations of the how-to manual approach. So you need to
    think of it as a case study and help someone decide the implications
    of their survey design.

Results:

-   see comment above

Discussion:

-   This starts just before “Research Opportunities”. You can keep these
    sections but there should be more preamble before jumping right into
    research opportunities. I suggest a general paragraph(s) that segue
    into your subsections of “research opportunities”, “future
    research”. In the Discussion you should address again some of the
    broader issues from the Introduction, e.g. how could spatial (depth)
    distribution changes anticipated under climate change for some
    population be tackled by survey design exploration now so that we
    can continue to track these changing populations 20 years down the
    road and do not lose the signal. How can this software help with
    that?

-   I think something that can be useful for readers is to outline the
    steps in a thought process a researcher might undertake when setting
    up a survey (perhaps a separate subsection) and then how one might
    go about a SimSurvey run for this. It also gives you a good
    opportunity to discuss your multispecies ideas:
    -   What is your current problem/needs
    -   What are the constraints (HR and $ perspective)
    -   What are the constraints from a biological, physical perspective
    -   Resources available
    -   What are your anticipated future needs (i.e. if climate change
        is going to make cod go deeper, will you be able to track it in
        20 years?)
    -   How could you consider some or all of these with SimSurvey
-   Your “Assumptions” section should come down into the Discussion

-   “Future Directions” should say something about randomfields re:
    Reviewer 1. Even if you just outline that you have considered it.
    You might also try to say something about optimisation of design
    which Reviewer 2 mentioned.
