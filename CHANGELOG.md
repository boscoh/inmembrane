# inmembrane changelog

## v0.94 (27-Sep-2013)
* Fixed LipoP and TMHMM web scrapers due to server side changes.
* TMB-HUNT plugin deprecated - server appears to be permanently offline
* Added result_poll_retries option to signalp_web, so test fails fast if job is stuck in PENDING state
* Made 'requests' v2.0.0 the lowest version for new installs
* Made default inmembrane.config use SOAP web services (*_web), included commented _scrape_web options in file.

## v0.93.2 (03-Dec-2012)
* Bugfixes to tmhmm_scrape_web and lipop_scrape_web plugins
* Handle Python 2.6 when OrderedDict is absent.

## v0.93 (03-Dec-2012)

* Added --version flag
* Report PSE loop length in gram_pos protocol
* Added tmhmm_scrape_web and lipop_scrape_web plugins (as default since CBS SOAP service seems currently broken)
* Added safe_seqid to be used by programs that munge or fail on seqids with funky punctuation

## v0.92 (22-Oct-2012)

* Started CHANGELOG.
