using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DE
{
	public class Random : System.Random
	{
		override public int Next(int max)
		{
			return (int)(NextDouble() * (double)max);
		}

		public double Next(double min, double max)
		{
			return min + NextDouble() * (max - min);
		}
	}
}
