# Missing references in `references.bib`

`main (1).tex` cites 20 keys; `references.bib` currently has 5. These 11 are
cited but not in the bib file, so the paper will show `[?]`/undefined-
citation warnings until they're added. Ready-to-paste BibTeX for each is
below (pulled from the last known-complete version of the file, in git
history at commit `ae7b328`).

## Fix this one first — key typo, not actually missing

Your current bib has a `ClaridaGaliGertler200` entry (note: no trailing
`0`), but the paper cites `\citep{ClaridaGaliGertler2000}`. Either rename
your existing entry's key to add the missing `0`, or paste the block below
as a new entry — don't do both, or you'll get a duplicate-key warning.

```bibtex
@article{ClaridaGaliGertler2000,
  author  = {Clarida, Richard and Gal{\'i}, Jordi and Gertler, Mark},
  title   = {Monetary policy rules and macroeconomic stability: Evidence and some theory},
  journal = {Quarterly Journal of Economics},
  year    = {2000}, volume = {115}, number = {1}, pages = {147--180}
}
```

## Genuinely missing (10 entries)

```bibtex
@article{dieboldyilmaz2012,
  author  = {Diebold, Francis X. and Yilmaz, Kamil},
  title   = {Better to give than to receive: Predictive directional measurement of volatility spillovers},
  journal = {International Journal of Forecasting},
  year    = {2012}, volume = {28}, number = {1}, pages = {57--66}
}

@article{ando2022,
  author  = {Ando, Tomohiro and Greenwood-Nimmo, Matthew and Shin, Yongcheol},
  title   = {Quantile connectedness: Modeling tail behavior in the topology of financial networks},
  journal = {Management Science}, year = {2022}, volume = {68}, number = {4}, pages = {2401--2431}
}

@article{balli2023,
  author  = {Balli, Faruk and Balli, Hatice O. and Dang, Tuan Hai Nam and Gabauer, David},
  title   = {Contemporaneous and lagged {R2} decomposed connectedness approach: New evidence from the energy futures market},
  journal = {Finance Research Letters}, year = {2023}, volume = {57}, pages = {104168}
}

@article{choishin2022,
  author  = {Choi, Ji-Eun and Shin, Dong Wan},
  title   = {Quantile correlation coefficient: a new tail dependence measure},
  journal = {Statistical Papers}, year = {2022}, volume = {63}, number = {4}, pages = {1075--1104}
}

@article{genizi1993,
  author  = {Genizi, Abraham},
  title   = {Decomposition of {R2} in multiple regression with correlated regressors},
  journal = {Statistica Sinica}, year = {1993}, volume = {3}, number = {2}, pages = {407--420}
}

@misc{shahzadetal,
  author  = {Shahzad, Syed Jawad Hussain and Ryan, Michael and Gabauer, David},
  title   = {Pseudo-quantile {R2} connectedness [companion methodology paper]},
  year    = {2026}, note = {Working paper}
}

@article{stenfors2026decomposing,
  title   = {Decomposing the rate of inflation: Forecast-based connectedness among CPI components},
  author  = {Stenfors, Alexis and Shabani, Mimoza and Gabauer, David and Toporowski, Jan},
  journal = {Economic Modelling},
  volume  = {163}, pages = {107708}, year = {2026}, publisher = {Elsevier}
}

@article{antonakakis2020,
  author  = {Antonakakis, Nikolaos and Chatziantoniou, Ioannis and Gabauer, David},
  title   = {Refined Measures of Dynamic Connectedness Based on Time-Varying Parameter Vector Autoregressions},
  journal = {Journal of Risk and Financial Management},
  year    = {2020}, volume = {13}, number = {4}, pages = {84}, publisher = {MDPI}
}

@misc{Powell2021,
  author = {Powell, Jerome H.},
  title  = {Monetary Policy in the Time of {COVID-19}},
  year   = {2021},
  note   = {Speech at ``Macroeconomic Policy in an Uneven Economy,'' economic policy
            symposium sponsored by the Federal Reserve Bank of Kansas City, Jackson Hole,
            Wyoming, August 27}
}

@misc{Powell2022JH,
  author = {Powell, Jerome H.},
  title  = {Monetary Policy and Price Stability},
  year   = {2022},
  note   = {Speech at ``Reassessing Constraints on the Economy and Policy,'' economic
            policy symposium sponsored by the Federal Reserve Bank of Kansas City, Jackson
            Hole, Wyoming, August 26}
}
```

## Not cited by `main (1).tex`, in case you want them back too

These were in the pre-trim `references.bib` but aren't currently cited
anywhere in the paper, so they weren't included above: `stockwatson2016`,
`boivin2009`, `pasten2020`, `ciccarelli2010`, `koenkermachado1999`,
`bonaccolto2025`, `taylor1993`, `orphanides2001`, `wuxia2016` (note:
`wu2016measuring`, same Wu–Xia paper under a different key, is already in
your current bib and *is* cited).

Full old version, for reference: `git show ae7b328:paper/references.bib`.

Delete this file once you've reconciled the bibliography.
