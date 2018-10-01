import pstats
ps = pstats.Stats('profile.txt')
ps.sort_stats('cumulative').print_stats(70)
