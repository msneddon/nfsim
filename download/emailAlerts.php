
<?php


function sendEmailMessage($to, $subject, $body)
{
	require_once "PEAR.php";
	require_once "Mail.php";
	
	$from = "NFsim Administrator <nfsim.super.user@gmail.com>";
	
	$host = "ssl://smtp.gmail.com";
	//$port = "587";
	$port = "465"; // this change seems to fix SMTP authentication issues with google servers -- srinivas
	$username = "nfsim.super.user@gmail.com";
	$password = "opennf4343";
	
	$headers = array ('From' => $from,
	  'To' => $to,
	  'Subject' => $subject);
	$smtp = Mail::factory('smtp',
	  array ('host' => $host,
	    'port' => $port,
	    'auth' => true,
	    'username' => $username,
	    'password' => $password));
	
	$mail = $smtp->send($to, $headers, $body);
	
	if (PEAR::isError($mail)) {
		  echo("<p> <font color=#FF0000>Error sending email! :" . $mail->getMessage() . "</font></p>");
	} //else {
	  //echo("<p>Message successfully sent!</p>");
	  //}
}

?>