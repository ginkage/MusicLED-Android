package com.ginkage.musicled;

import android.os.SystemClock;
import android.util.Log;
import java.io.*;
import java.net.InetAddress;
import java.net.Socket;

public class TCPClient {
	public static final String SERVERIP = "10.10.100.254"; //your computer IP address
	public static final int SERVERPORT = 8899;
	private OnMessageReceived mMessageListener = null;
	private boolean mRun = false;

	private BufferedOutputStream out = null;
	private BufferedInputStream in = null;

	/**
	 *  Constructor of the class. OnMessagedReceived listens for the messages received from server
	 */
	public TCPClient(OnMessageReceived listener) {
		mMessageListener = listener;
	}

	/**
	 * Sends the message entered by client to the server
	 * @param msg text entered by client
	 */
	public void sendMessage(int[] msg) {
		if (out != null) {
			try {
				for (int i = 0; i < msg.length; i++)
					out.write(msg[i]);
				out.flush();
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void stopClient() {
		mRun = false;
	}

	public boolean isConnected() {
		return (in != null);
	}

	public void run() {
		mRun = true;

		try {
			//here you must put your computer's IP address.
			InetAddress serverAddr = InetAddress.getByName(SERVERIP);

			Log.i("TCP Client", "C: Connecting...");

			//create a socket to make the connection with the server
			Socket socket = new Socket(serverAddr, SERVERPORT);

			try {
				//send the message to the server
				out = new BufferedOutputStream(socket.getOutputStream());

				//receive the message which the server sends back
				in = new BufferedInputStream(socket.getInputStream());

				Log.i("TCP Client", "C: Connected.");

				if (mMessageListener != null)
					mMessageListener.messageReceived("Connected.");

				//in this while the client listens for the messages sent by the server
				int[] msg = new int[] { 0xFF };
				while (mRun) {
					sendMessage(msg);
					SystemClock.sleep(1000);
				}
			}
			catch (Exception e) {
				Log.e("TCP", "S: Error", e);
				in = null;
				out = null;
				mRun = false;
				if (mMessageListener != null)
					mMessageListener.messageReceived(e.getMessage());
			}
			finally {
				//the socket must be closed. It is not possible to reconnect to this socket
				// after it is closed, which means a new socket instance has to be created.
				socket.close();
				Log.i("TCP Client", "C: Disonnected.");
			}
		}
		catch (Exception e) {
			Log.e("TCP", "C: Error", e);
			in = null;
			out = null;
			mRun = false;
			if (mMessageListener != null)
				mMessageListener.messageReceived(e.getMessage());
		}
	}

	//Declare the interface. The method messageReceived(String message) will must be implemented in the MyActivity
	//class at on asynckTask doInBackground
	public interface OnMessageReceived {
		public void messageReceived(String message);
	}
}
