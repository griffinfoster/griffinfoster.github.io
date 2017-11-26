---
layout: post
title:  "ALFABURST: A commensal search for Fast Radio Bursts with Arecibo"
date:   2017-11-26
excerpt: "ALFABURST: A commensal search for Fast Radio Bursts with Arecibo"
tag:
- python
- FRBs
- Arecibo
comments: false
---

In November 2017 the ALFABURST group published our [two-year results](https://arxiv.org/abs/1710.10806) in our survey to detect now fast radio bursts (FRBs) using the Arecibo radio telescope. This is a survey that runs commensally with all observations using the ALFA multi-beam feed. This has allowed for a search for weak (possibly very distant) FRB sources which is not possible with less sensitive telescopes. In that time we did not detect any new FRB sources.

<img src="https://griffinfoster.github.io/assets/img/ALFABURST/sensitivity_range.png" height="300">

*ALFABURST single pulse sensitivity (purple region). Previously detected FRBs from Parkes (black triangle), GBT (red circle), Arecibo (white diamond), UTMOST (teal pentagon), and ASKAP (yellow-green hexagon) are plotted for reference. Line of constant fluence (solid) are plotted for reference. The fluence completeness (dashed) is 0.5 Jy ms out to pulse widths of 16 ms.*

Presently, bright FRBs are being detected by wide field of view telescopes. But, there have been a limited number of detections of weaker FRBs. This is in part due to the more sensitive telescopes having a smaller field of view. As there is so far no preferential direction to find new FRBs, the most likely way to detect new FRBs is to look in as much sky as possible, i.e. a large survey speed. The limited detection of weak FRBs could also be due to the a change in the number of events. The predicted number of FRBs is based on the dispersion measure-distance relation.

<img src="https://griffinfoster.github.io/assets/img/ALFABURST/cartview_sky_coverage.png" height="400">

*ALFABURST survey coverage, most of the survey coverage is during the PALFA and AGES survey, shown in boxes.*

The main work of our published manuscript is that we have processed two years worth of detections by our automated detection pipeline. Over 250k events were recorded, the vast majority of which were due to radio frequency interference (RFI) and the system switching between observing modes. To process all the events, a sub-sample was labelled by hand into 9 classes of events, mainly different types of RFI and systematics, and pulses from previously detected pulsars. A random-forest prioritizer model was trained on the labelled datasets using approximately 400 extracted features. This model created a queue of events to examine based on the likelihood that the detected event was due to a pulse.

An unknown pulse (shown below) was detected while the telescope was moving between fields. This is likely a pulse from an unknown pulsar. Follow-up observations are in progress.

<img src="https://griffinfoster.github.io/assets/img/ALFABURST/Beam5_fb_D20170618T005616_buffer2_spectrum.png" height="400">

With this non-detection result, we have been able to put limit on the number of weak FRBs. The survey will continue to operate while ALFA is in use, and we will be upgrading our instrumentation to process more bandwidth, increasing our sensitivity.