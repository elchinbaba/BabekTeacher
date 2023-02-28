using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DE
{
	public class Ranges
	{
		public readonly double[] RMin;
		public readonly double[] RMax;

		public Ranges(double[] rmin, double[] rmax)
		{
			if (rmin != null && rmax != null)
			{
				RMin = rmin;
				RMax = rmax;
			}
			else
			{
				new Ranges(0.0, 1.0);
			}
		}

		public Ranges(double rmin, double rmax)
		{
			RMin = new double[] { rmin };
			RMax = new double[] { rmax };
		}

	}
}
