using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DE
{
	public class SearchControls
	{
		protected internal int PopSize;
		protected internal double F;
		protected internal double Cr;

		public SearchControls(int ps, double f, double cr)
		{
			PopSize = ps;
			F = f;
			Cr = cr;
			if (PopSize < 10) PopSize = 10;
		}

		public SearchControls(int ps)
		{
			PopSize = ps;
			F = 0.5;
			Cr = 0.9;
			if (PopSize < 10) PopSize = 10;
		}

		public SearchControls()
		{
			PopSize = 60;
			F = 0.5;
			Cr = 0.9;
		}

	}
}
