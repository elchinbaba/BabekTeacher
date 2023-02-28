namespace DE
{
	public class BestStrategy : SearchStrategy
	{
		override public void Apply(
			double f, double cr, int n, double[] x, double[] best,
			double[] r1, double[] r2, double[] r3
		)
		{
			int j = RNG.Next(n);
			for (int counter = 0; counter < n; counter++)
			{
				if (RNG.NextDouble() < cr) x[j] = best[j] + f * (r2[j] - r3[j]);
				j = (j + 1) % n;
			}
		}
	}
}
