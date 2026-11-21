# Does the timing result survive a change in pulse energy?

At the reference source setting, sending both pulses early allows charged daughters to appear. Moving one part in 1,600 of the emitted energy to a later time keeps the channel closed. This comparison persists to wall time 11 on all three grids. Changing the original pulse dose by just 0.1% in either direction also closes the channel in the all-early trial through time 11. Thus the statement about one part in 1,600 describes the reference pulse. A larger total dose does not guarantee particle production in this model. [Numerical results](CHECK.json)

| Source dose relative to reference | Total emitted energy on original grid | Peak scalar, both early | Peak scalar, one late | Which case produced daughters by time 11? |
|---|---:|---:|---:|---|
| 99.9% | 0.97705948 | 0.358825 | 0.356282 | Neither |
| 100% | 0.97803744 | 0.748040 | 0.360361 | Both early only |
| 100.1% | 0.97901540 | 0.527008 | 0.527205 | Neither |

The decay channel opens at scalar value 0.53. The middle all-early trial produced its first daughters near 6.256. By wall time 11 their cumulative rest mass was 0.0021203 in model units. Both quantities are internal simulation outputs. The two columns of peak values include evolution after the first daughters form whenever the channel opens.

Within each row we matched the emitted energy by calibrating the delayed pulse against the evolving outer surface. Over all nine matched pairs the largest relative energy difference was below $2.5\times10^{-13}$ on the reference machine. Each pair was run on the original grid, at half the original time step, and at twice the original radial cell counts. All reached the common wall time and agreed on which cases produced daughters. The peak scalar for the successful case changed from 0.74804 to 0.74425 on doubling the radial resolution. The unsuccessful higher-dose, all-early case peaked between 0.527007 and 0.527009 across grids.

Only three doses were tested. They bracket a successful setting without measuring the width or the number of successful intervals. The model has several stages in which fields exchange energy. The higher-dose trials reach the final release region and turn back before opening. Distinguishing a shift in timing from a change in scalar work in that region requires a finer dose scan and a time-resolved comparison of the individual energy transfers. The outer pulse shape and angularly averaged gravitational model are inherited from the [source calculation](../source/README.md).

To repeat this comparison, point the script at the unpacked parent and at the output of the source calculation.

~~~sh
python -I -B sensitivity/energy.py --parent /path/to/replayed-parent \
    --source /path/to/SOURCE.json --output /tmp/pulse-dose.json --jobs 2
~~~
