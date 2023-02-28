using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DE
{
	public class AdaptiveControls : SearchControls
	{
		protected Random RNG;

		virtual public void Update()
		{
			//
			F = (RNG.Next(0.0, 1.0) < 0.1) ? RNG.Next(0.0, 1.0) : 0.5;
			Cr = (RNG.Next(0.0, 1.0) < 0.1) ? RNG.Next(0.5, 1.0) : 0.9;
			/*
			F= 0.9;
			Cr= 1.0;
			*/
		}

		public AdaptiveControls() : base()
		{
			RNG = new Random();
		}

		public AdaptiveControls(Problem dep) : this()
		{
			PopSize = dep.NPars * 10;
		}

		public AdaptiveControls(int ps) : this()
		{
			PopSize = ps;
		}

		public AdaptiveControls(Optimizer deo) : this()
		{
			PopSize = deo.DEP.NPars * 10;
		}

	}
}
