namespace ZFIS
{
	public class _ZNumber3
	{
		public readonly FuzzySet A;
		public readonly FuzzySet B;

		public const double Infinity = (+1.0) / (+0.0);
		public const double Neginfinity = (+1.0) / (-0.0);

		public _ZNumber3(double v1, double v2, double v3, double v4, double v5, double v6)
		{
			if (Min(v4, v5, v6) < 0 || Max(v4, v5, v6) > 1)
			{
				System.Console.WriteLine("Check Z-number's B-part format");
				System.Console.ReadKey();
				throw new System.FormatException();
			}
			//A= new FuzzySet(new double[]{ v1, v2, v3 }, new double[]{ 0.00000001, 1, 0.00000001 });
			//B= new FuzzySet(new double[]{ v4, v5, v6 }, new double[]{ 0.00000001, 1, 0.00000001 });
			A = new FuzzySet(
				new double[] { v1, 0.8 * v1 + 0.2 * v2, 0.6 * v1 + 0.4 * v2, 0.4 * v1 + 0.6 * v2, 0.2 * v1 + 0.8 * v2, v2, 0.8 * v2 + 0.2 * v3, 0.6 * v2 + 0.4 * v3, 0.4 * v2 + 0.6 * v3, 0.2 * v2 + 0.8 * v3, v3 },
				new double[] { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0 }
			);
			B = new FuzzySet(
				new double[] { v4, 0.8 * v4 + 0.2 * v5, 0.6 * v4 + 0.4 * v5, 0.4 * v4 + 0.6 * v5, 0.2 * v4 + 0.8 * v5, v5, 0.8 * v5 + 0.2 * v6, 0.6 * v5 + 0.4 * v6, 0.4 * v5 + 0.6 * v6, 0.2 * v5 + 0.8 * v6, v6 },
				new double[] { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0 }
			);
		}

		/*
		public double AL{ get{ return A.Value[0]; } }
		public double AM{ get{ return A.Value[1]; } }
		public double AR{ get{ return A.Value[2]; } }
		public double BL{ get{ return B.Value[0]; } }
		public double BM{ get{ return B.Value[1]; } }
		public double BR{ get{ return B.Value[2]; } }
		*/
		public double AL { get { return A.Value[0]; } }
		public double AM { get { return A.Value[5]; } }
		public double AR { get { return A.Value[10]; } }
		public double BL { get { return B.Value[0]; } }
		public double BM { get { return B.Value[5]; } }
		public double BR { get { return B.Value[10]; } }

		public _ZNumber3(double v1, double v2, double v3) :
		this(v1, v2, v3, 1, 1, 1)
		{
		}

		public _ZNumber3(double v1, double v2, double v3, double v4) :
		this(v1, (v1 + v2) / 2, v2, v3, (v3 + v4) / 2, v4)
		{
		}

		public _ZNumber3(double v1, double v2) :
		this(v1, (v1 + v2) / 2, v2, 1, 1, 1)
		{
		}

		public _ZNumber3(double v) :
		this(v, v, v, 1, 1, 1)
		{
		}

		double Max(double x1, double x2, double x3)
		{
			if (x1 >= x2) { if (x1 >= x3) return x1; }
			else { if (x2 >= x3) return x2; }
			return x3;
		}

		double Max(double x1, double x2)
		{
			return x1 >= x2 ? x1 : x2;
		}

		double Min(double x1, double x2)
		{
			return x1 <= x2 ? x1 : x2;
		}

		double Min(double x1, double x2, double x3)
		{
			if (x1 <= x2) { if (x1 <= x3) return x1; }
			else { if (x2 <= x3) return x2; }
			return x3;
		}

		public _Interval Membership(double x)
		{
			if (x >= AL && x <= AM && AM - AL < 1E-8 || x >= AM && x <= AR && AR - AM < 1E-8) return new _Interval(BL, BR);
			double amemb = System.Math.Max(0, System.Math.Min((x - AL) / (AM - AL), (AR - x) / (AR - AM)));
			return new _Interval(amemb * BL, amemb * BR);
		}

		override public string ToString()
		{
			return "[ " + AL + ", " + AM + ", " + AR + " ] [ " + BL + ", " + BM + ", " + BR + " ]";
		}

		public _Interval ToInterval()
		{
			return new _Interval(AL, AR);
		}

		public double ToDouble()
		{
			return ToInterval().ToDouble();
		}

		public FuzzySet ToFuzzySet()
		{
			return A;
		}
	}
}
