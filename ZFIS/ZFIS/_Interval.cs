namespace ZFIS
{
	public class _Interval
	{
		public readonly double L;
		public readonly double R;

		public _Interval(double l, double r)
		{
			L = System.Math.Min(l, r);
			R = System.Math.Max(l, r);
		}

		public _Interval(double d)
		{
			L = d;
			R = d;
		}

		public static bool operator <(_Interval a, _Interval b)
		{
			return (a.L + a.R < b.L + b.R) ? true : false;
		}

		public static bool operator >(_Interval a, _Interval b)
		{
			return (a.L + a.R > b.L + b.R) ? true : false;
		}

		public static bool operator <=(_Interval a, _Interval b)
		{
			return (a.L + a.R <= b.L + b.R) ? true : false;
		}

		public static bool operator >=(_Interval a, _Interval b)
		{
			return (a.L + a.R >= b.L + b.R) ? true : false;
		}

		public double ToDouble()
		{
			return (L + R) / 2;
		}

		override public string ToString()
		{
			return "[ " + L + ", " + R + " ]";
		}
	}
}
