
This folder is a sort of haphazard collection of scripts that we used in the
paper. Not all of them ended up providing results that went directly into the
paper, and a lot of them were also absorbed into the extensive testing in the
package itself. But we include them all here for the curious people who would
like to obtain results from the paper. Of course, as the versions of BesselK.jl
change the exact numbers here will also change a bit, although if I did
everything right the first tagged release of this code (or maybe the initial
commit) should be almost the exact source that was used to generate the results
in v1 of the paper.

Should anything come up that you'd like to discuss, please don't hesitate to
contact me. You can find my email addresses on my website, which you can find by
googling my name (Chris Geoga).

Some misc notes:

-- I would _not_ suggest using the fitting scripts in `./demo/` as the basis of
your own code for estimating parameters. You could certainly do much worse and
it does leverage my generic go-to package `GPMaxlik.jl`, which has a ton of nice
features (not that I'm biased). But I'd sooner suggest looking at the example
files in that repo as a template.

-- I've re-organized the code a bit and some code lives in `../examples/`. I've
done my best to make sure that all of these tests still run as they should after
the re-org, but if something doesn't, please open an issue or email me or
something. It's probably just an accident that I will be able to resolve
immediately.

