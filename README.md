# Observing Strategy Task Force

#### Charter
To respond to the white paper call (https://www.lsst.org/call-whitepaper-2018) on behalf of the DESC.

----

## Deliverables
-Two white papers, one focused on WFD and one focused on DDF and other mini-Surveys, to be led by Michelle Lochner and Dan Scolnic.  These white papers are due Nov. 30.  The overleaf for DDF is here: https://www.overleaf.com/read/xgpxqdmzydcc
(publically viewable) and the overleaf for WFD is here: https://www.overleaf.com/read/ywbstqjtnctk.  

-One journal article that details all the metrics established for the white papers.  Expectation for timeline to be roughly one month behind the white papers.  

-While not a set deliverable, there is an expectation that individual contributions will form their own papers which can go into significant detail of the analysis that relates to observing strategy.


----------------------------------


## Calendar

July 19, Task Force Meeting.  Presentation of additional information needed for call response.

July 24, CMU ObsStrat XWG Session at 11:00.  July 26, CMU ObsStrat XWG Session at 09:00.

July 25, Settled list of additional information needed for appropriate response to call

August 2, Task Force Meeting

August 10, First assessment of subset of WFD and mini-survey strategies, with metrics.  Slides from each WG.

August 13-16, LSST Community Workshop including multiple ObsStrat sessions

August 16, Task Force Meeting

August 30, Task Force Meeting

September 17-19, LSST Cadence Hackathon at Flatiron (https://www.lsstcorporation.org/node/193).  Need to apply.

September 21, Full assessment of strategies and metrics due for task force

Sep 27 - moved to weekly meetings every Thursday at 9am Central.

Oct 11 - First drafts of white papers shared with Task Force.

* [ ] **October 25**, Revised white papers sent out to DESC reviewers.

* [ ] **October 30**, DESC Seminar presentation on White Papers.

* [ ] **November 5**, Responses from DESC due.

* [ ] **November 30, !!!!White Papers Due!!!**

The key takeaways are a first take due August 10 before the community workshop so we can present to the Project any issues we might have, and a full assessment by September 15.  Than iterations over the next two months with white papers sent out to DESC reviewers on Oct 25.


----------------------------------

## White Paper Template


The LSST Project wants to standardize responses to its call. The latex template provided is here: https://github.com/lsst-pst/survey_strategy_wp/blob/master/WP_submission_template.tex#L84.
Its sections are as follows:

> a. Abstract —
> Please provide a short summary of your scientific goals and survey strategy modifications here.
> 
> b. General Information —
> Author Info, Science Category, Survey Type Category, Observing Strategy Category
> 
> c.  Scientific Motivation  —
> Describe the scientific justification for this white paper in the context of your field, as well as the importance to the general program of astronomy, including the relevance over the next decade. Describe other relevant data, and justify why LSST is the best facility for these observations.
>
> d.  Technical Description —
> Describe your survey strategy modifications or proposed observations. Please comment on each observing constraint below, including the technical motivation behind any constraints. Where relevant, indicate if the constraint applies to all requested observations or a specific subset. Please note which constraints are not relevant or important for your science goals.
> 
> e. Performance Evaluation -
> Please describe how to evaluate the performance of a given survey in achieving your desired science goals, ideally as a heuristic tied directly to the observing strategy (e.g. number of visits obtained within a window of time with a specified set of filters) with a clear link to the resulting effect on science.
>
> f. Special Data Processing --
> Describe any data processing requirements beyond the standard LSST Data Management pipelines and how these will be achieved.

[Updated after talking to Lynne Jones] Here is what we feel we need to add to this template, and what we will be asking of each WG:

Add:  A section on the strategies considered, including detailed information about all the information we needed to add to the white paper call (e.g., locations for different fields) that would be assessed by each group.

Add from each group a 2-3 page description:

  i.  Brief introduction of motivation of each probe and what is unique to LSST
  
  ii. How analysis was done/ code was used [Please link to Jupyter notebook]
  
  iii.  Justification of the metric
 
  iv.  Explanation of preferences of working group with qualifiers about how good/bad strategies are

  v.  What is needed to know if people outside group use the metric as is without consultation

  vi.  How would external data change rankings/metric?
  
  vii. What is further work / slight modifications to survey strategy thought about?
  
Add: Summary of how external (non-LSST) can be used, and what changes in terms of final rankings/preferences assuming external data


----

## Summary of New OpSim Runs
The details of the new OpSim runs are here: https://docushare.lsst.org/docushare/dsweb/Get/Document-28716 

Below they are listed as a handy summary (thanks Philippe!):

1. baseline2018a -  Project-official baseline (official opsim v4 baseline, 6/2018). No dome crawl.
2. kraken_2026  - Unofficial baseline (expected next baseline, including dome crawl). 
3. kraken_2035 - Survey with 9 DD mini surveys instead of 5. Addition of 4 extra extragalactic DD fields. 
4. kraken_2036 - A rolling cadence in WFD. Three regions of declination, alernating on/off every third year from years 3-8. Standard WFD coverage in years 1-3, 8-10 (first two and last two years). 
5. colossus_2665 - Baseline style, with slightly expanded WFD footprint to better satisfy SRD constraints on fO 
6. colossus_2664 - No separated Galactic Plane proposal. Simply run WFD over galactic plane region. 
7. colossus_2667 - Survey with single visits only. Standard visits in pairs each night replaced with single visits per night. 
8. pontus_2002 - Very large WFD footprint (24,700 sq deg) and 5 DD mini-surveys (5 fields).
9. pontus_2489 - Survey with more visits. Standard 2x15s visits replaced with 1x20s in grizy and 1x40s in u band. 
10. pontus_2502 - A rolling cadence in WFD. Two regions of declination, alternating on/off every other year throughout the entire survey, but the background WFD remains "on" at a 25% reward level. 
11. mothra_2045 - A rolling cadence in WFD. Two regions of declination, alternating on/off every other year throughout entire survey. 


## Additional Information
In the above, we mention additional information.  By that we mean that there are a number of pieces of information we as a group need to fill in about the call.  We are asking the following people to prepare this for our next meeting:

1.  Location of WFD - Humna Awan.

2.  Location of DDFs - David Rubin.

3.  List of WFD strategies not in white paper call - Michelle Lochner

4.  Cadence of DDF - Philippe Gris

5.  Special calibration observations - Nicolas Regnault

6.  ToO Capability - Marcelle Soares-Santos

1-4 will likely go in section C above.  5 and 6 may have their own dedicated place within the DDF and mini-surveys response. 

For 1-4, we are asking these people to come up with a contained menu of options that each group can consider as part of their response.  We would like to have this settled by the CMU meeting, with presentation on July 19.

-----

## Working with other surveys
* For each of our white papers, we will have a dedicated section on the benefits of external datasets.  This will include ones from Euclid, WFIRST and other spectroscopic surveys.  
* Our process will be totally transparent and our github page will be public - there is nothing that needs to be hidden about our deliberations.  However, we want to maintain a ‘unified front’ within DESC, so we ask for everyone to refrain from stating any conclusions of the task force until these conclusions have been reached.
* Endorsement of strategies outside of DESC will be done within the DESC white papers.

## DESC work:
- documentation of the data products such as dither strategies and simlibs [created for DESC work](doc/README.md)
