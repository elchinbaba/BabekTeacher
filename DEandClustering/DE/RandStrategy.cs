using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DE
{
	public class RandStrategy : SearchStrategy
	{
		override public void Apply(
			double f, double cr, int n, double[] x, double[] best,
			double[] r1, double[] r2, double[] r3
		)
		{
			int j = RNG.Next(n);
			for (int counter = 0; counter < n; counter++)
			{
				if (RNG.NextDouble() < cr) x[j] = r1[j] + f * (r2[j] - r3[j]);
				j = (j + 1) % n;
			}
		}
	}
}
