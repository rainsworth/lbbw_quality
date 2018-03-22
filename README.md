# lbbw_quality

Quality Assessment working group for the LOFAR Long Baseline Pipeline during the Lorentz Centre workshop on [High Resolution Surveying with International LOFAR](https://www.lorentzcenter.nl/lc/web/2018/983/info.php3?wsid=983&venue=Snellius), from 19 Mar 2018 through 23 Mar 2018.

We are adapting quality control scripts (astrometry, flux accuracy, dynamic range, etc.) from the [DDFacet](https://github.com/mhardcastle/ddf-pipeline) [Quality Pipeline](https://github.com/mhardcastle/ddf-pipeline/blob/master/scripts/quality_pipeline.py) (which has a [GNU General Public License](https://github.com/mhardcastle/ddf-pipeline/blob/master/LICENSE.md)) for use in the LOFAR [Long Baseline Pipeline](https://github.com/lmorabit/long_baseline_pipeline) (which doesn't yet have a License).


We are currently testing on CEP3. To run:

```
$ module load lofar
$ module load prefactor
$ genericpipeline.py lbbw_quality/quality_lb.parset -d -c pipeline.cfg
```

(pipeline.cfg is required for the generic pipeline so is included here for testing - it should be removed when integrating this repo into the Long Baseline Pipeline as it will already have this file.)


### Contributors

* [Alexander Kappes](https://github.com/alexmatze)
* [Frits Sweijen](https://github.com/tikk3r)
* Maria Arias
* [Rachael Ainsworth](https://github.com/rainsworth)
* [Sean Mooney](https://github.com/mooneyse)
