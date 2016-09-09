


<?php


function verifyEmail ($tryEmail, $tryPassword, &$valid)
{
	//we are always a no go until confirmation
	$valid=0;
		
	//connect to DB
	// old code: mysql_connect('localhost','nfWebRegUser','nepo4nEFsim') or die('Error 1.');
	// what follows is new code for the new web site -- srinivas
	mysql_connect('mysql.emonet.biology.yale.edu','nfwebreguser','nepo4nEFsim') or die('Error 1.');
	//@mysql_select_db('NFsimWebRegister') or die('Error 2.');
	//new database is all lowercase
	@mysql_select_db('nfsimwebregister2') or die('Error 2.');

	//get data
	$getUserQuery="SELECT passwd FROM nfRegUser where email=\"$tryEmail\";";
	$names=mysql_query($getUserQuery);
	$rows=mysql_numrows($names);
	
	//search data (could be made faster)
	$i=0;
	while($i<$rows) {
		$passwd = mysql_result($names,$i,"passwd");
		
		if(strcmp($passwd,$tryPassword) == 0) {
			$valid=1;
			//echo("found<br>");
			break;	
		}
		
		$i++;	
	};
	
	
	//exit database
	mysql_close();
	
	//return
	return;
} 

?>