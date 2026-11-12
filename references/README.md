# Published RE1 coordinate regression evidence

`re1_coordinate_validation.json` contains five full human reference oligos from
[Mouri et al., Table S2](https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001926-mmc3.zip),
joined by explicit oligo identifier to the authors'
[FIMO output](https://github.com/tewhey-lab/duoREST/blob/main/processed_data/ref_all_fimo.txt).
Their [analysis](https://github.com/tewhey-lab/duoREST/blob/main/README.md#L180)
selects MA0138.2 hits at **one-based inclusive 91–111**. The corresponding Python
slice is `[90:111]`, and this repository's zero-based inclusive span is `(90, 110)`.

The examples include both motif orientations, strong sites, a weak site, and a
site immediately below the 20.86 threshold. Full oligos, matched motif bases,
orientation, and published scores provide independent regression evidence; the
test does not merely repeat a configured constant. `make pilot-validate` checks
these examples with the production matrix and scorer before publishing its
comparison table. The installed-FIMO integration test uses the same matrix
copied to `tests/fixtures/MA0138.2.meme` (JASPAR MA0138.2, MD5
`4de4483fff043233a0f077d2731d7275`), so it does not require research downloads.

These reference records test coordinate convention and scoring compatibility.
They do not establish binding, evolutionary direction, or gate behavior.
