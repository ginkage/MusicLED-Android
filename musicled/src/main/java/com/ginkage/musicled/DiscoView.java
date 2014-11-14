package com.ginkage.musicled;

import android.content.Context;
import android.graphics.*;
import android.media.audiofx.Visualizer;
import android.os.AsyncTask;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.util.AttributeSet;
import android.view.SurfaceHolder;
import android.view.SurfaceView;
import android.view.View;
import android.widget.TextView;

/**
 * View that draws, takes keystrokes, etc. for a simple DiscoLander game.
 *
 * Has a mode which RUNNING, PAUSED, etc. Has a x, y, dx, dy, ... capturing the
 * current ship physics. All x/y etc. are measured with (0,0) at the lower left.
 * updatePhysics() advances the physics based on realtime. draw() renders the
 * ship, and does an invalidate() to prompt another draw() as soon as possible
 * by the system.
 */
public class DiscoView extends SurfaceView implements SurfaceHolder.Callback {
	class DiscoThread extends Thread {
		/*
		 * Member (state) fields
		 */
		/**
		 * Current height of the surface/canvas.
		 *
		 * @see #setSurfaceSize
		 */
		private int mCanvasHeight = 1;

		/**
		 * Current width of the surface/canvas.
		 *
		 * @see #setSurfaceSize
		 */
		private int mCanvasWidth = 1;

		/** Message handler used by thread to interact with TextView */
		private Handler mHandler;

		/** Indicate whether the surface has been created & is ready to draw */
		private boolean mRun = false;

		private final Object mRunLock = new Object();

		/** Handle to the surface manager object we interact with */
		private SurfaceHolder mSurfaceHolder;

		public float fps = 0;
		private long start_frame;
		private long frames_drawn;

		private Visualizer mVisualizer;
		private int mCaptureSize;
		private int mSamplingRate;
		private byte[] fft;
		private byte[] wave = new byte[65536];
		private double[] left = new double[16384];
		private double[] right = new double[16384];

		private static final int winSize = 6;
		private double[] winR = new double[winSize];
		private double[] winG = new double[winSize];
		private double[] winB = new double[winSize];
		private int winIdx = 0;

		private boolean hifi = false;

		private final int M = 14;
		private final int N = (1 << M);
		private FFT fftCalc = new FFT(M);
		private double[] re = new double[N];
		private double[] im = new double[N];

        private int drawMode = 0;

		public DiscoThread(SurfaceHolder surfaceHolder, Context context, Handler handler)
		{
			// get handles to some important objects
			mSurfaceHolder = surfaceHolder;
			mHandler = handler;
			mContext = context;

			mVisualizer = new Visualizer(0);
			int[] sizeRange = mVisualizer.getCaptureSizeRange();
			mVisualizer.setCaptureSize(sizeRange[1]);
			int res = mVisualizer.setScalingMode(2); //Visualizer.SCALING_MODE_NORMALIZED);

			mCaptureSize = mVisualizer.getCaptureSize();
			mSamplingRate = mVisualizer.getSamplingRate();
			fft = new byte[mCaptureSize];

			if (res == Visualizer.SUCCESS) {
				mCaptureSize *= 16;
				hifi = true;
			}
		}

		@Override
		public void run() {
			while (mRun) {
				Canvas c = null;
				try {
					c = mSurfaceHolder.lockCanvas(null);
					synchronized (mSurfaceHolder) {
						// Critical section. Do not allow mRun to be set false until
						// we are sure all canvas draw operations are complete.
						//
						// If mRun has been toggled false, inhibit canvas operations.
						synchronized (mRunLock) {
							if (mRun && c != null) doDraw(c);
						}
					}
				} finally {
					// do this in a finally so that if an exception is thrown
					// during the above, we don't leave the Surface in an
					// inconsistent state
					if (c != null) {
						mSurfaceHolder.unlockCanvasAndPost(c);
					}
				}
			}
		}

		/**
		 * Used to signal the thread whether it should be running or not.
		 * Passing true allows the thread to run; passing false will shut it
		 * down if it's already running. Calling start() after this was most
		 * recently called with false will result in an immediate shutdown.
		 *
		 * @param b true to run, false to shut down
		 */
		public void setRunning(boolean b) {
			// Do not allow mRun to be modified while any canvas operations
			// are potentially in-flight. See doDraw().
			synchronized (mRunLock) {
				mRun = b;
				mVisualizer.setEnabled(b);
			}
		}

		public void setMessage(CharSequence message) {
			synchronized (mSurfaceHolder) {
				Message msg = mHandler.obtainMessage();
				Bundle b = new Bundle();
				b.putString("text", message.toString());
				b.putInt("viz", View.VISIBLE);
				msg.setData(b);
				mHandler.sendMessage(msg);
			}
		}

		/* Callback invoked when the surface dimensions change. */
		public void setSurfaceSize(int width, int height) {
			// synchronized to make sure these all change atomically
			synchronized (mSurfaceHolder) {
				mCanvasWidth = width;
				mCanvasHeight = height;
				start_frame = System.currentTimeMillis();
			}
		}

		private double clamp(double val)
		{
			if (val < 0) val += 12;
			if (val > 6) val = 12 - val;
			return val / 6.0;
		}

        public void toggleDisplay()
        {
            drawMode = (drawMode + 1) % 3;
        }

		private void doDraw(Canvas canvas)
		{
			long curTime = System.currentTimeMillis();
			if (curTime > start_frame + 1000) {
				fps = frames_drawn * 1000.0f / (curTime - start_frame);
				start_frame = curTime;
				frames_drawn = 0;
				setMessage(String.format("%.1f FPS", fps));
			}

			if (hifi) {
                mVisualizer.getWaveForm(wave);

                int k;
                for (k = 0; k < mCaptureSize; k++) {
                    short l = (short)((wave[4*k + 1] << 8) | wave[4*k + 0]);
                    short r = (short)((wave[4*k + 3] << 8) | wave[4*k + 2]);
                    left[k] = (l / 32768.0);
                    right[k] = (r / 32768.0);
                }

                if (drawMode == 0)
				    doDrawHiFiWave(canvas);
                else
                    doDrawHiFiFFT(canvas);
            }
			else
				doDrawWave(canvas);

			frames_drawn++;
		}

		private void doDrawHiFiFFT(Canvas canvas)
		{
			int k;
			for (k = 0; k < mCaptureSize; k++) {
				re[k] = (left[k] + right[k]) * 0.5;
				im[k] = 0;
			}

			fftCalc.fft(re, im);

            if (drawMode == 1)
			    canvas.drawRGB(0, 0, 0);
			Paint p = new Paint();

			int maxFreq = mCaptureSize;
			double base = Math.log(Math.pow(2, 1.0 / 12.0));

			double fcoef = Math.pow(2, 57.0 / 12.0) / 440.0; // Frequency 440 is a note number 57 = 12 * 4 + 9
			double minFreq = mSamplingRate / maxFreq;
			double minNote = Math.log(minFreq * fcoef) / base;
			double minOctave = Math.floor(minNote / 12.0);
			fcoef = Math.pow(2, ((4 - minOctave) * 12 + 9) / 12.0) / 440.0; // Shift everything by several octaves
			double maxNote = Math.log(mSamplingRate * fcoef) / base;

			int baseY = (mCanvasHeight * 3) / 4;

			double kx = mCanvasWidth / maxNote;
			double ky = mCanvasHeight * 0.5;
			int lastx = -1;
			int maxY = 0;
			
			double maxAmp = 0;
			int maxR = 0, maxG = 0, maxB = 0;

			for (k = 1; k < maxFreq >> 1; k++) {
				double amp = Math.hypot(re[k], im[k]) / 256.0;

                if (drawMode == 1) {
                    double frequency = (k * (double) mSamplingRate) / maxFreq;
                    double note = Math.log(frequency * fcoef) / base; // note = 12 * Octave + Note

                    int x = (int) Math.round(note * kx);
                    int y = (int) Math.round(amp * ky);

                    if (y > maxY)
                        maxY = y;

                    if (lastx != x) {
                        double spectre = note - 12.0 * Math.floor(note / 12.0); // spectre is within [0, 12)

                        lastx = x;
                        maxY = 0;

                        double R = clamp(spectre - 6);
                        double G = clamp(spectre - 10);
                        double B = clamp(spectre - 2);

                        double mx = Math.max(Math.max(R, G), B);
                        double mn = Math.min(Math.min(R, G), B);
                        double mm = mx - mn;
                        if (mm == 0) mm = 1;

                        R = (R - mn) / mm;
                        G = (G - mn) / mm;
                        B = (B - mn) / mm;

                        p.setARGB(255, (int) Math.round(R * 255), (int) Math.round(G * 255), (int) Math.round(B * 255));
                        canvas.drawLine(x, baseY, x, baseY - y, p);
                    }
                }
				else if (amp > maxAmp) {
					maxAmp = amp;
                    double frequency = (k * (double) mSamplingRate) / maxFreq;
                    double note = Math.log(frequency * fcoef) / base; // note = 12 * Octave + Note
					double spectre = note - 12.0 * Math.floor(note / 12.0); // spectre is within [0, 12)

	                int ring = (96 - (int) Math.floor(spectre * 8)); // [1 .. 96]
	                int cring = ring + 48;
	                if (cring > 100)
		                cring -= 100; // [1 .. 100]
	                setColor(cring);

					double R = clamp(spectre - 6);
					double G = clamp(spectre - 10);
					double B = clamp(spectre - 2);

					double mx = Math.max(Math.max(R, G), B);
					double mn = Math.min(Math.min(R, G), B);
					double mm = mx - mn;
					if (mm == 0) mm = 1;

					maxR = (int) Math.round(255.0 * (R - mn) / mm);
					maxG = (int) Math.round(255.0 * (G - mn) / mm);
					maxB = (int) Math.round(255.0 * (B - mn) / mm);
				}
			}
            if (drawMode != 1)
			    canvas.drawRGB(maxR, maxG, maxB);
		}

		private void doDrawHiFiWave(Canvas canvas)
		{
			canvas.drawRGB(0, 0, 0);

			Paint p = new Paint();
			p.setARGB(255, 255, 255, 255);

            int k;
			int lastx = 0;
			int minL = 0, maxL = 0, minR = 0, maxR = 0;
			int baseY = mCanvasHeight / 4;

			for (k = 0; k < mCaptureSize; k++) {
				int x = (int) Math.round((mCanvasWidth * k) / (double)(mCaptureSize - 1));
				int ry = (int) Math.round(left[k] * mCanvasHeight * 0.25);
				int ly = (int) Math.round(right[k] * mCanvasHeight * 0.25);

				if (ly < minL) minL = ly;
				if (ly > maxL) maxL = ly;
				if (ry < minR) minR = ry;
				if (ry > maxR) maxR = ry;

				if (x != lastx) {
					lastx = x;
					canvas.drawLine(x, baseY - minL, x, baseY - maxL, p);
					canvas.drawLine(x, baseY*3 - minR, x, baseY*3 - maxR, p);
					minL = maxL = minR = maxR = 0;
				}
			}
		}

		private void doDrawWave(Canvas canvas)
		{
			mVisualizer.getWaveForm(wave);
			Paint p = new Paint();
			p.setARGB(255, 255, 255, 255);
			canvas.drawRGB(0, 0, 0);

			int baseY = mCanvasHeight / 2;
			int k;

			for (k = 0; k < mCaptureSize; k++) {
				int d = (wave[k] & 0xFF);
				double h = (d - 128) / 255.0;

				int x = (int) Math.round((mCanvasWidth * k) / (double)(mCaptureSize - 1));
				int y = (int) Math.round(baseY + h * mCanvasHeight * 0.5);

				canvas.drawLine(x, baseY, x, y, p);
			}
		}

		private void doDrawFFT(Canvas canvas)
		{
			mVisualizer.getFft(fft);

			Paint p = new Paint();

			int k, maxFreq = mCaptureSize / 2;
			double base = Math.log(Math.pow(2, 1.0 / 12.0));

			double fcoef = Math.pow(2, 57.0 / 12.0) / 440.0; // Frequency 440 is a note number 57 = 12 * 4 + 9
			double minFreq = mSamplingRate / maxFreq;
			double minNote = Math.log(minFreq * fcoef) / base;
			double minOctave = Math.floor(minNote / 12.0);
			fcoef = Math.pow(2, ((4 - minOctave) * 12 + 9) / 12.0) / 440.0; // Shift everything by several octaves
			double maxNote = Math.log(mSamplingRate * fcoef) / base;

			double maxPeak = 0;
			int maxR = 0, maxG = 0, maxB = 0;
			double sumR = 0, sumG = 0, sumB = 0;

			canvas.drawRGB(0, 0, 0);

			int baseY = (mCanvasHeight * 3) / 4;

			for (k = 1; k < maxFreq; k++) {
				double frequency = (k * (double) mSamplingRate) / maxFreq;
				double real = fft[k*2];
				double imag = fft[k*2 + 1];
				double amp = Math.hypot(real, imag) / 255.0;

				double note = Math.log(frequency * fcoef) / base; // note = 12 * Octave + Note
				double spectre = note - 12.0 * Math.floor(note / 12.0); // spectre is within [0, 12)

				double R = clamp(spectre - 6);
				double G = clamp(spectre - 10);
				double B = clamp(spectre - 2);

				double mx = Math.max(Math.max(R, G), B);
				double mn = Math.min(Math.min(R, G), B);
				double mm = mx - mn;
				if (mm == 0) mm = 1;

				R = (R - mn) / mm;
				G = (G - mn) / mm;
				B = (B - mn) / mm;

				sumR += R * amp;
				sumG += G * amp;
				sumB += B * amp;

				int x = (int) Math.round((mCanvasWidth * note) / maxNote); //* (double) (k - 1)) / (maxFreq - 1));
				int y = (int) Math.round(amp * mCanvasHeight * 0.5);

				p.setARGB(255, (int) Math.round(R * 255), (int) Math.round(G * 255), (int) Math.round(B * 255));
				canvas.drawLine(x, baseY, x, baseY - y, p);

				if (amp > maxPeak)
					maxPeak = amp;
		  }
/*
			winR[winIdx] = sumR;
			winG[winIdx] = sumG;
			winB[winIdx] = sumB;

			int i;
			sumR = sumG = sumB = 0;
			for (i = 0; i < winSize; i++) {
				sumR += winR[i];
				sumG += winG[i];
				sumB += winB[i];
			}

			double mx = Math.max(Math.max(sumR, sumG), sumB);
			double mn = Math.min(Math.min(sumR, sumG), sumB);
			double mm = mx - mn;
			if (mm == 0) mm = 1;

			int R = (int) Math.round(255.0 * maxPeak * (sumR - mn) / mm);
			int G = (int) Math.round(255.0 * maxPeak * (sumG - mn) / mm);
			int B = (int) Math.round(255.0 * maxPeak * (sumB - mn) / mm);

			canvas.drawRGB(R, G, B);
*/
			winIdx = (winIdx + 1) % winSize;
		}
	}

	/** Handle to the application context, used to e.g. fetch Drawables. */
	private Context mContext;

	/** Pointer to the text view to display "Paused.." etc. */
	private TextView mStatusText;

	/** The thread that actually draws the animation */
	private DiscoThread thread;

	private TCPClient mTcpClient = null;

	public DiscoView(Context context, AttributeSet attrs) {
		super(context, attrs);

		// register our interest in hearing about changes to our surface
		SurfaceHolder holder = getHolder();
		holder.addCallback(this);

		mContext = context;

		setFocusable(true); // make sure we get key events
	}

	public class connectTask extends AsyncTask<String,String,TCPClient> {
		@Override
		protected TCPClient doInBackground(String... message) {
			//we create a TCPClient object and
			mTcpClient = new TCPClient(new TCPClient.OnMessageReceived() {
				@Override
				//here the messageReceived method is implemented
				public void messageReceived(String message) {
					//this method calls the onProgressUpdate
					publishProgress(message);
				}
			});
			mTcpClient.run();

			return null;
		}

		@Override
		protected void onProgressUpdate(String... values) {
			super.onProgressUpdate(values);
		}
	}

	private void openSocket() {
		if (mTcpClient == null || !mTcpClient.isConnected())
			new connectTask().execute("");
	}

	private void closeSocket() {
		if (mTcpClient != null) {
			mTcpClient.stopClient();
			mTcpClient = null;
		}
	}

	private void sendPacket(int[] msg) {
		if (mTcpClient != null && mTcpClient.isConnected())
			mTcpClient.sendMessage(msg);
	}

	private void sendMessage(int[] msg) {
		if (msg.length == 1)
			sendPacket(msg);
		else if (msg.length == 4) {
			int[] newMsg = new int[] { 0x55, 0x34, 0x33, 0x39, 0x02, 0x00, 0, 0, 0, 0, 0xAA, 0xAA };
			newMsg[6] = msg[0];
			newMsg[7] = msg[1];
			newMsg[8] = msg[2];
			newMsg[9] = msg[3];
			sendPacket(newMsg);
		}
	}

	private void setColor(int color)
	{
		color -= 4;
		if (color >= 1 && color <= 96) {
			int[] msg = new int[] { 0x01, 0x01, color, color + 4 };
			sendMessage(msg);
		}
	}

	/**
	 * Fetches the animation thread corresponding to this DiscoView.
	 *
	 * @return the animation thread
	 */
	public DiscoThread getThread() {
		return thread;
	}

	/**
	 * Installs a pointer to the text view used for messages.
	 */
	public void setTextView(TextView textView) {
		mStatusText = textView;
	}

	/* Callback invoked when the surface dimensions change. */
	public void surfaceChanged(SurfaceHolder holder, int format, int width, int height) {
		thread.setSurfaceSize(width, height);
	}

	/*
	 * Callback invoked when the Surface has been created and is ready to be
	 * used.
	 */
	public void surfaceCreated(SurfaceHolder holder) {
		// start the thread here so that we don't busy-wait in run()
		// waiting for the surface to be created

		// create thread only; it's started in surfaceCreated()
		thread = new DiscoThread(holder, mContext, new Handler() {
			@Override
			public void handleMessage(Message m) {
				//noinspection ResourceType
				mStatusText.setVisibility(m.getData().getInt("viz"));
				mStatusText.setText(m.getData().getString("text"));
			}
		});

		thread.setRunning(true);
		thread.start();
		thread.setMessage("Dummy Message");
		openSocket();
	}

	/*
	 * Callback invoked when the Surface has been destroyed and must no longer
	 * be touched. WARNING: after this method returns, the Surface/Canvas must
	 * never be touched again!
	 */
	public void surfaceDestroyed(SurfaceHolder holder) {
		// we have to tell thread to shut down & wait for it to finish, or else
		// it might touch the Surface after we return and explode
		boolean retry = true;
		thread.setRunning(false);
		while (retry) {
			try {
				thread.join();
				retry = false;
			} catch (InterruptedException e) {
			}
		}

		thread = null;
		closeSocket();
	}

    public void toggleDisplay()
    {
        if (thread != null)
            thread.toggleDisplay();
    }

	public class FFT
	{
		int n, m;

		// Lookup tables.  Only need to recompute when size of FFT changes.
		double[] cos;
		double[] sin;
		double[] window;

		public FFT(int m)
		{
			int n = 1 << m;

			this.n = n;
			this.m = m;

			// precompute tables
			cos = new double[n / 2];
			sin = new double[n / 2];

			for (int i = 0; i < n / 2; i++) {
				cos[i] = Math.cos(-2 * Math.PI * i / n);
				sin[i] = Math.sin(-2 * Math.PI * i / n);
			}

			// Make a blackman window:
			window = new double[n];
			for (int i = 0; i < window.length; i++)
				window[i] = 0.42 - 0.5 * Math.cos(2 * Math.PI * i / (n - 1)) + 0.08 * Math.cos(4*Math.PI*i/(n-1));
		}

		/***************************************************************
		 * fft.c
		 * Douglas L. Jones
		 * University of Illinois at Urbana-Champaign
		 * January 19, 1992
		 * http://cnx.rice.edu/content/m12016/latest/
		 *
		 *   fft: in-place radix-2 DIT DFT of a complex input
		 *
		 *   input:
		 * n: length of FFT: must be a power of two
		 * m: n = 2**m
		 *   input/output
		 * x: double array of length n with real part of data
		 * y: double array of length n with imag part of data
		 *
		 *   Permission to copy and use this program is granted
		 *   as long as this header is included.
		 ****************************************************************/
		public void fft(double[] x, double[] y)
		{
			int i, j, k, n1, n2, a;
			double c, s, e, t1, t2;


			// Bit-reverse
			j = 0;
			n2 = n / 2;
			for (i = 1; i < n - 1; i++) {
				n1 = n2;
				while (j >= n1) {
					j = j - n1;
					n1 = n1 / 2;
				}

				j = j + n1;

				if (i < j) {
					t1 = x[i];
					x[i] = x[j];
					x[j] = t1;
					t1 = y[i];
					y[i] = y[j];
					y[j] = t1;
				}
			}

			// FFT
			n1 = 0;
			n2 = 1;

			for (i = 0; i < m; i++) {
				n1 = n2;
				n2 = n2 + n2;
				a = 0;

				for (j=0; j < n1; j++) {
					c = cos[a];
					s = sin[a];
					a += 1 << (m - i - 1);

					for (k = j; k < n; k = k + n2) {
						t1 = c * x[k + n1] - s * y[k + n1];
						t2 = s * x[k + n1] + c * y[k + n1];
						x[k + n1] = x[k] - t1;
						y[k + n1] = y[k] - t2;
						x[k] = x[k] + t1;
						y[k] = y[k] + t2;
					}
				}
			}
		}
	}
}
