function mailme(subject, message)
if nargin == 1
    message = '';
end

destination = char('pnfs/ib{po' - 1);
send_text_message(destination, 'gmail', subject, message);